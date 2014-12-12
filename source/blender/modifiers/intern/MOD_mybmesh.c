/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/modifiers/intern/MOD_mybmesh.c
 *  \ingroup modifiers
 */

/* This code is based of the tessellation part of the paper 
 * "Computing Smooth Surface Contours with Accurate Topology"
 * (Pierre Bénard, Aaron Hertzmann, Michael Kass).
 * Currently available at:
 * http://www.labri.fr/perso/pbenard/publications/contours.html
 *
 * The numbers in the comments refers to the chapters in the paper.
 */


#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"

#include "BLI_math.h"
#include "BLI_string.h"
#include "BLI_utildefines.h"
#include "BLI_buffer.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_modifier.h"
#include "BKE_deform.h"
#include "BKE_subsurf.h"

#include "bmesh.h"
#include "bmesh_tools.h"

#include "MOD_util.h"

#include "intern/CCGSubSurf.h"

//TODO this modifier depends on OSD. So if it's not compiled in, remove this modifier
#include <opensubdiv/osdutil/evaluator_capi.h>

struct OpenSubdiv_EvaluatorDescr;

typedef struct {
	BMEdge *orig_edge;
	BMFace *orig_face;
	float u, v;
} Vert_buf;

typedef struct {
	BMEdge *cusp_e;
	float cusp_co[3];
} Cusp;

//TODO for Kr look in subdiv.cpp in coutours source code (II)

//TODO dynamic arrays, use BLI_stack, BLI_buffer, BLI_mempool, BLI_memarena.

static void verts_to_limit(BMesh *bm, struct OpenSubdiv_EvaluatorDescr *eval){
	
	int i, j;
	
	BMIter iter_v, iter_f;
	BMVert *vert;
    BMFace *f;

    //TODO is it possible to only get non adjusted verts?
	//IE not moving a vert more than once.

	BM_ITER_MESH_INDEX (f, &iter_f, bm, BM_FACES_OF_MESH, i) {
			BM_ITER_ELEM_INDEX (vert, &iter_v, f, BM_VERTS_OF_FACE, j) {
				float u,v, du[3], dv[3];
				switch(j){
					case 1 : u = 1, v = 0;
							 break;
					case 2 : u = v = 1;
							 break;
					case 3 : u = 0, v = 1;
							 break;
					default: u = v = 0;
							 break;
				}
				openSubdiv_evaluateLimit(eval, i, u, v, vert->co, du, dv);
				//Adjust vert normal to the limit normal
				cross_v3_v3v3(vert->no, dv, du);
				normalize_v3(vert->no);
				//printf("j: %d\n",j);
			}
			//printf("i: %d\n",i);
			//printf("face i: %d\n", BM_elem_index_get(f));
	}

}

static bool calc_if_B(const float cam_loc[3], const float P[3], const float du[3], const float dv[3]){
	//Is the point back facing?
	float nor[3], view_vec[3];

	cross_v3_v3v3(nor, dv, du);
	//TODO normalization is probably not needed
	normalize_v3(nor);
	sub_v3_v3v3(view_vec, cam_loc, P);

	return ( dot_v3v3(nor, view_vec) < 0);
}

static bool calc_if_B_nor(const float cam_loc[3], const float P[3], const float nor[3]){
	//Is the point back facing?
	float view_vec[3];

	sub_v3_v3v3(view_vec, cam_loc, P);

	return ( dot_v3v3(nor, view_vec) < 0);
}


static float get_facing_dir(const float cam_loc[3], const float P[3], const float du[3], const float dv[3]){
	//Get if point is back facing (-) or front facing (+). Zero if it's directly on the contour
	float nor[3], view_vec[3];

	cross_v3_v3v3(nor, dv, du);
	//TODO normalization is probably not needed
	normalize_v3(nor);
	sub_v3_v3v3(view_vec, cam_loc, P);

	return dot_v3v3(nor, view_vec);
}


static float get_facing_dir_nor(const float cam_loc[3], const float P[3], const float nor[3]){
	//Get if point is back facing (-) or front facing (+). Zero if it's directly on the contour
	float view_vec[3];

	sub_v3_v3v3(view_vec, cam_loc, P);

	return dot_v3v3(nor, view_vec);
}

static void split_edge_and_move_vert(BMesh *bm, BMEdge *edge, const float new_pos[3],
									const float du[3], const float dv[3], const float u, const float v){
	//Split edge one time and move the created vert to new_pos
	
    BMVert *vert;
	printf("Split edge!\n");

	BM_mesh_elem_hflag_disable_all(bm, BM_EDGE, BM_ELEM_TAG, false);
	BM_elem_flag_enable(edge, BM_ELEM_TAG);
	BMO_op_callf(bm, BMO_FLAG_DEFAULTS,
			"subdivide_edges edges=%he cuts=%i quad_corner_type=%i use_single_edge=%b",
			BM_ELEM_TAG, 1, SUBD_STRAIGHT_CUT, true);
	//Get the newly created vertex
	vert = BM_vert_at_index_find(bm, BM_mesh_elem_count(bm, BM_VERT) - 1);
	copy_v3_v3(vert->co, new_pos);

	//Adjust vert normal to the limit normal
	cross_v3_v3v3(vert->no, dv, du);
	normalize_v3(vert->no);

	//TODO save what face the new vert belonged to. Also u,v coords
}


static void split_edge_and_move_cusp(BMesh *bm, BMEdge *edge, const float new_pos[3]){
	//Split edge one time and move the created vert to new_pos
	
    BMVert *vert;
	printf("Cusp split!\n");

	BM_mesh_elem_hflag_disable_all(bm, BM_EDGE, BM_ELEM_TAG, false);
	BM_elem_flag_enable(edge, BM_ELEM_TAG);
	BMO_op_callf(bm, BMO_FLAG_DEFAULTS,
			"subdivide_edges edges=%he cuts=%i quad_corner_type=%i use_single_edge=%b",
			BM_ELEM_TAG, 1, SUBD_STRAIGHT_CUT, true);
	//Get the newly created vertex
	vert = BM_vert_at_index_find(bm, BM_mesh_elem_count(bm, BM_VERT) - 1);
	copy_v3_v3(vert->co, new_pos);
}
	
static void get_uv_coord(BMVert *vert, BMFace *f, float *u, float *v){
	//Get U,V coords of a vertex
    int i;
	BMIter iter;
	BMVert *vert_iter;

	BM_ITER_ELEM_INDEX (vert_iter, &iter, f, BM_VERTS_OF_FACE, i) {
		if(vert == vert_iter){
			switch(i){
				case 1 : *u = 1, *v = 0;
						 break;
				case 2 : *u = *v = 1;
						 break;
				case 3 : *u = 0, *v = 1;
						 break;
				default: *u = *v = 0;
						 break;
			}
		}
	}
}

static void split_BB_FF_edges(BMesh *bm, BMesh *bm_orig, struct OpenSubdiv_EvaluatorDescr *eval,
							  const float cam_loc[3], BLI_Buffer *new_vert_buffer){
    //Split BB,FF edges if they have sign crossings
	
	int i, face_index;
	BMIter iter_f, iter_e;
	BMEdge *e;
	BMFace *f;
	BMVert *v1, *v2;
	float v1_u, v1_v, v2_u, v2_v, step;
	bool is_B;
	int orig_edges = BM_mesh_elem_count(bm_orig, BM_EDGE);
    int initial_edges = BM_mesh_elem_count(bm, BM_EDGE);

	//Do 10 samples but don't check end and start point
	step = 1.0f/12.0f;

	BM_ITER_MESH_INDEX (e, &iter_e, bm, BM_EDGES_OF_MESH, i) {
        Vert_buf v_buf;

        if( !(i < initial_edges) ){
			//Now we are working on edges we added in this function
			break;
		}

        //TODO perhaps we need to use the limit surface normal
		
        is_B = calc_if_B_nor(cam_loc, e->v1->co, e->v1->no);

        if( is_B  != calc_if_B_nor(cam_loc, e->v2->co, e->v2->no) ){
			//This is not a FF or BB edge
			continue;
		}

		if( i < orig_edges ){
            //This edge exists on the original mesh
			//TODO why do I have to use find? Segfault otherwise...
			//remember to replace the rest of "at_index"
			//why is table_ensure not fixing the assert?
			BMEdge *orig_e = BM_edge_at_index_find(bm_orig, i);

			//Get face connected to edge from orig mesh
			//TODO is it wise to use BM_ITER_ELEM here?
			BM_ITER_ELEM (f, &iter_f, orig_e, BM_FACES_OF_EDGE) {
                //Get first face
				break;
			}     
			
			face_index = BM_elem_index_get(f);

			v1 = orig_e->v1;
			v2 = orig_e->v2;

			v_buf.orig_edge = orig_e;
			v_buf.orig_face = NULL;
		} else {
			BMVert *vert_arr[2];
			bool found_face;

			//This should be safe because the vert count is still the same as the original mesh.
			v1 = BM_vert_at_index_find( bm_orig, BM_elem_index_get( e->v1 ) ); 
			v2 = BM_vert_at_index_find( bm_orig, BM_elem_index_get( e->v2 ) ); 

            vert_arr[0] = v1;
			vert_arr[1] = v2;

            //TODO add checks if to hande if there is no face
			found_face = BM_face_exists_overlap(vert_arr, 2, &f);
			face_index = BM_elem_index_get(f);

			v_buf.orig_edge = NULL;
			v_buf.orig_face = f;
		}

		//TODO can I just check each vert once?
		get_uv_coord(v1, f, &v1_u, &v1_v);
		get_uv_coord(v2, f, &v2_u, &v2_v);

		//TODO perhaps we need to check the limit normal and not the vertex normal?
		/*
		{
			float P1[3], P2[3], du[3], dv[3];
			//is this a FF or BB edge?
			openSubdiv_evaluateLimit(eval, face_index, v1_u, v1_v, P1, du, dv);

			is_B = calc_if_B(cam_loc, P1, du, dv);

			openSubdiv_evaluateLimit(eval, face_index, v2_u, v2_v, P2, du, dv);

			if( is_B  != calc_if_B(cam_loc, P2, du, dv) ){
				//FB edge, we only want to split FF or BB
				//Skip to next edge
				continue;  
			}
		}
        */
		{
            int i;
            float u, v;
			float P[3], du[3], dv[3];

			if( v1_u == v2_u ){
				u = v1_u;
				v = step;

				for(i=0; i < 10; i++){
					openSubdiv_evaluateLimit(eval, face_index, u, v, P, du, dv);
					if( calc_if_B(cam_loc, P, du, dv) != is_B ){
						split_edge_and_move_vert(bm, e, P, du, dv, u, v);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(new_vert_buffer, Vert_buf, v_buf);
						break;
					}
					v += step;
				}
			} else if ( v1_v == v2_v ){
				u = step;
				v = v1_v;

				for(i=0; i < 10; i++){
					openSubdiv_evaluateLimit(eval, face_index, u, v, P, du, dv);
					if( calc_if_B(cam_loc, P, du, dv) != is_B ){
						split_edge_and_move_vert(bm, e, P, du, dv, u, v);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(new_vert_buffer, Vert_buf, v_buf);
						break;
					}
					u += step;
				}
			} else {
				float step_u;
				if((v1_u == 0 && v1_v == 0) || (v2_u == 0 && v2_v == 0)){
					step_u = step;
					u = v = step;
				} else {
					step_u = -step;
					u = 1.0f - step;
					v = step;
				}
				for(i=0; i < 10; i++){
					openSubdiv_evaluateLimit(eval, face_index, u, v, P, du, dv);
					if( calc_if_B(cam_loc, P, du, dv) != is_B ){
						split_edge_and_move_vert(bm, e, P, du, dv, u, v);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(new_vert_buffer, Vert_buf, v_buf);
						break;
					}
					u += step_u;
					v += step;
				}
			}
		}

	}

}

static float get_k_r(struct OpenSubdiv_EvaluatorDescr *eval, int face_index, float u, float v, const float cam_loc[3]){
	float du[3], dv[3], P[3], no[3], d1[3], d2[3];
    float k1, k2;
	float I[2][2], II[2][2];
	openSubdiv_evaluateLimit(eval, face_index, u, v, P, du, dv);

	cross_v3_v3v3(no, dv, du);
	normalize_v3(no);

    //http://en.wikipedia.org/wiki/Principal_curvature

	I[0][0] = dot_v3v3(du, du);
	I[0][1] = dot_v3v3(du, dv);
	I[1][0] = dot_v3v3(du, dv);
	I[1][1] = dot_v3v3(dv, dv);

	//get dudu dudv dvdv
	{
		float dudu[3], dudv[3], dvdv[3];
		float du_old[3], dv_old[3];

        float step_u, step_v;

		copy_v3_v3(du_old, du);
		copy_v3_v3(dv_old, dv);

		if( u < 0.9f ){
			step_u = 0.1f;
		} else {
			step_u = -0.1f;
		}

		openSubdiv_evaluateLimit(eval, face_index, u + step_u, v, P, du, dv);

        sub_v3_v3v3(dudu, du, du_old);
        sub_v3_v3v3(dudv, dv, dv_old);

		mul_v3_fl(dudu, 1.0f/step_u);
		mul_v3_fl(dudv, 1.0f/step_u);

		if( v < 0.9f ){
			step_v = 0.1f;
		} else {
			step_v = -0.1f;
		}

		openSubdiv_evaluateLimit(eval, face_index, u, v + step_v, P, du, dv);

        sub_v3_v3v3(dvdv, dv, dv_old);
		mul_v3_fl(dvdv, 1.0f/step_v);

		II[0][0] = dot_v3v3(dudu, no);
		II[0][1] = dot_v3v3(dudv, no);
		II[1][0] = dot_v3v3(dudv, no);
		II[1][1] = dot_v3v3(dvdv, no);

		copy_v3_v3(du, du_old);
		copy_v3_v3(dv, dv_old);

	}
	{
		float S[2][2];
		float detI = determinant_m2(I[0][0], I[0][1], I[1][0], I[1][1]);

		if(fabsf(detI) < 1e-20){
			detI = 1e-20;
			printf("detI near zero!!!\n");
		}

		S[0][0] = (II[0][1]*I[0][1] - II[0][0]*I[1][1]) / detI;
		S[0][1] = (II[1][1]*I[0][1] - II[0][1]*I[1][1]) / detI;
		S[1][0] = (II[0][0]*I[0][1] - II[0][1]*I[0][0]) / detI;
		S[1][1] = (II[0][1]*I[0][1] - II[1][1]*I[0][0]) / detI;

		{
			//TODO perhaps remove pdir2
			float pdir1[2], pdir2[2];
			float traceS = S[0][0] + S[1][1];
			float detS = determinant_m2(S[0][0], S[0][1], S[1][0], S[1][1]);
			float diff = traceS*traceS - 4.0f * detS;

			if(diff >= 0){
				float sqrtDiff = sqrt3f(diff);
				k1 = 0.5f * (traceS + sqrtDiff);
				k2 = 0.5f * (traceS - sqrtDiff);
				if(fabsf(k1) < fabsf(k2)){
					float swap = k1;
					k1 = k2;
					k2 = swap;
				}
				if(fabsf(S[1][0]) > 1e-20){
					copy_v2_fl2(pdir1, k1 - S[1][1], S[1][0]);
					copy_v2_fl2(pdir2, k2 - S[1][1], S[1][0]);
				}else if (fabsf(S[0][1]) > 1e-20){
					copy_v2_fl2(pdir1, S[0][1], k1 - S[0][0]);
					copy_v2_fl2(pdir2, S[0][1], k2 - S[0][0]);
				}
				normalize_v2(pdir1);
				normalize_v2(pdir2);
			}else{
				//k1 = 0.0f;
				//k2 = 0.0f;
				//d1 = tanU;
				//d2 = tanV;
				printf("diff neg\n");
				return 0;
			}

			mul_v3_fl(du, pdir1[0]);
			mul_v3_fl(dv, pdir1[1]);
            add_v3_v3v3(d1, du, dv);
            normalize_v3(d1);
			cross_v3_v3v3(d2, no, d1);
            //d2 = d1 ^ limitNormal; //tanU * pdir2[0] + tanV * pdir2[1];
		}
	}
	{
		float view_vec[3], ndotv, sintheta, u2, v2, k_r;

		sub_v3_v3v3(view_vec, P, cam_loc);
		normalize_v3(view_vec);

		ndotv = dot_v3v3(no, view_vec);
		sintheta = 1.0f - ndotv*ndotv;
		u = dot_v3v3(view_vec, d1);
		v = dot_v3v3(view_vec, d2);
		u2 = u*u;
		v2 = v*v;
		k_r = (k1 * u2 + k2 * v2) / sintheta;
		return k_r;
	}
}

static void convert_uv_to_new_face(BMEdge *e, BMFace *f, float *u, float *v){
	//convert the old u/v coords to the new face coord
	
	float v1_u, v1_v, v2_u, v2_v, step;
	float old_u = *u;
	get_uv_coord(e->v1, f, &v1_u, &v1_v);
	get_uv_coord(e->v2, f, &v2_u, &v2_v);

	//where did we split this edge? (uv coords)
	if( old_u == 0 || old_u == 1){
		step = *v;
	} else {
		step = *u;
	}

	if( v1_u == v2_u ){
		*u = v1_u;
        *v = step;
	} else if( v1_v == v2_v ){
        *u = step;
		*v = v1_v;
	} else {
		//This should not happen...
		printf("convert_uv_to\n");
	}
	
}

static void bisect_search(const float v1_uv[2], const float v2_uv[2], struct OpenSubdiv_EvaluatorDescr *eval,
						BMesh *bm, BMEdge *e, int face_index, const float cam_loc[3]){
	//Search edge for sign crossing and split it!
	int i;
	float face_dir, uv_P[2], P[3], du[3], dv[3];
	float step = 0.5f;
	float step_len = 0.25f;
	float v1_face = get_facing_dir_nor(cam_loc, e->v1->co, e->v1->no);
	for( i = 0; i < 10; i++){
		interp_v2_v2v2( uv_P, v1_uv, v2_uv, step);		
		openSubdiv_evaluateLimit(eval, face_index, uv_P[0], uv_P[1], P, du, dv);
        face_dir = get_facing_dir(cam_loc, P, du, dv);

        if( face_dir == 0 ){
			//We got lucky and found the zero crossing!
			break;
		}

		if( (face_dir < 0) == (v1_face < 0) ){
			step += step_len;
		} else {
			step -= step_len;
		}
		step_len = step_len/2.0f;
	}
	
	if( len_v3v3(P, e->v1->co) < 1e-3 || len_v3v3(P, e->v2->co) < 1e-3 ){
		//Do not insert a new vert here
		return;
	}

	split_edge_and_move_vert(bm, e, P, du, dv, uv_P[0], uv_P[1]);
}

static void contour_insertion(BMesh *bm, BMesh *bm_orig, BLI_Buffer *new_vert_buffer,
							  BLI_Buffer *cusp_edges, struct OpenSubdiv_EvaluatorDescr *eval,
							  const float cam_loc[3]){
    int i, cusp_i, face_index;
	BMEdge *e;
	BMFace *f;
	BMIter iter_e;
	float v1_u, v1_v, v2_u, v2_v;

	int orig_verts = BM_mesh_elem_count(bm_orig, BM_VERT);
	int orig_edges = BM_mesh_elem_count(bm_orig, BM_EDGE);
    int initial_edges = BM_mesh_elem_count(bm, BM_EDGE);

    printf("Buffer count: %d\n", new_vert_buffer->count);

	BM_ITER_MESH_INDEX (e, &iter_e, bm, BM_EDGES_OF_MESH, i) {
        if( !(i < initial_edges) ){
			//Now we are working on edges we added in this function
			break;
		}

        //Check if this is a cusp edge
		for(cusp_i = 0; cusp_i < cusp_edges->count; cusp_i++){
			Cusp cusp = BLI_buffer_at(cusp_edges, Cusp, cusp_i);
			if(cusp.cusp_e == e){
				//Do not split this edge yet
				printf("skipped cusp edge\n");
				continue;
			}
		}

        //TODO perhaps we need to use the limit surface normal
        if( calc_if_B_nor(cam_loc, e->v1->co, e->v1->no) == calc_if_B_nor(cam_loc, e->v2->co, e->v2->no) ){
			//This is not a FB or BF edge
			continue;
		}

		{
			int v1_idx = BM_elem_index_get(e->v1); 
			int v2_idx = BM_elem_index_get(e->v2);
			Vert_buf v_buf1, v_buf2;
			BMVert *v1 = NULL, *v2 = NULL;
			bool v1_has_face = false, v2_has_face = false;

			if( (v1_idx + 1) > orig_verts){
				printf("num1: %d\n", v1_idx - orig_verts);
				v_buf1 = BLI_buffer_at(new_vert_buffer, Vert_buf, v1_idx - orig_verts);
				v1_u = v_buf1.u;
				v1_v = v_buf1.v;
				if( v_buf1.orig_face ){
					v1_has_face = true;
				}
			} else {
				v1 = BM_vert_at_index_find( bm_orig, v1_idx ); 
			}
			if( (v2_idx + 1) > orig_verts){
				printf("num2: %d\n", v2_idx - orig_verts);
				v_buf2 = BLI_buffer_at(new_vert_buffer, Vert_buf, v2_idx - orig_verts);
				v2_u = v_buf2.u;
				v2_v = v_buf2.v;
				if( v_buf2.orig_face ){
					v2_has_face = true;
				}
			} else {
				v2 = BM_vert_at_index_find( bm_orig, v2_idx ); 
			}

			if( v1 && v2 ){
				printf("v1 & v2\n");
				//TODO make this a external function
				if( i < orig_edges ){ 
				//this edge is on the original mesh
                BMIter iter_f;

				//TODO why do I have to use find? Segfault otherwise...
				//remember to replace the rest of "at_index"
				//why is table_ensure not fixing the assert?
				BMEdge *orig_e = BM_edge_at_index_find(bm_orig, i);

				//Get face connected to edge from orig mesh
				//TODO is it wise to use BM_ITER_ELEM here?
				BM_ITER_ELEM (f, &iter_f, orig_e, BM_FACES_OF_EDGE) {
					//Get first face
					break;
				}
				} else {
					BMVert *vert_arr[] = {v1 ,v2};
					BM_face_exists_overlap(vert_arr, 2, &f);
				}
				get_uv_coord(v1, f, &v1_u, &v1_v);
				get_uv_coord(v2, f, &v2_u, &v2_v);
			} else if ( v1 ){
				printf("v1\n");
				if( v2_has_face ){
					printf("face\n");
					f = v_buf2.orig_face; 
				} else {
					BMVert *vert_arr[3];

					vert_arr[0] = v_buf2.orig_edge->v1;
					vert_arr[1] = v_buf2.orig_edge->v2;
                    vert_arr[2] = v1;
                    //TODO check if get face fails
					if( v_buf2.orig_edge->v1 == v2 || v_buf2.orig_edge->v2 == v2 ){
						printf("Same vert!\n");
						BM_face_exists_overlap(vert_arr, 2, &f);
					} else {
						BM_face_exists_overlap(vert_arr, 3, &f);
					}
					convert_uv_to_new_face( v_buf2.orig_edge, f, &v2_u, &v2_v);
				}
				get_uv_coord(v1, f, &v1_u, &v1_v);
			} else if ( v2 ){
				printf("v2\n");
				if( v1_has_face ){
					printf("face\n");
					f = v_buf1.orig_face; 
				} else {
					BMVert *vert_arr[3];

					vert_arr[0] = v_buf1.orig_edge->v1;
					vert_arr[1] = v_buf1.orig_edge->v2;
                    vert_arr[2] = v2;
                    //TODO check if get face fails
					if( v_buf1.orig_edge->v1 == v2 || v_buf1.orig_edge->v2 == v2 ){
						printf("Same vert!\n");
						BM_face_exists_overlap(vert_arr, 2, &f);
					} else {
						BM_face_exists_overlap(vert_arr, 3, &f);
					}
					convert_uv_to_new_face( v_buf1.orig_edge, f, &v1_u, &v1_v);
				}
				get_uv_coord(v2, f, &v2_u, &v2_v);
			} else {
                printf("Two verts on the same orig face\n");
				if( v1_has_face || v2_has_face ){
					if( v1_has_face && v2_has_face ){
						printf("face 1 & 2\n");
						f = v_buf1.orig_face; 
					} else if ( v1_has_face ) {
						printf("face1\n");
						f = v_buf1.orig_face; 
						convert_uv_to_new_face( v_buf2.orig_edge, f, &v2_u, &v2_v);
					} else {
						printf("face2\n");
						f = v_buf2.orig_face; 
						convert_uv_to_new_face( v_buf1.orig_edge, f, &v1_u, &v2_v);
					}
				} else {
					//No orig face. So this in on a orig edge. So just get the face from the v1 edge
					BMVert *vert_arr[] = {v_buf1.orig_edge->v1 ,v_buf1.orig_edge->v2};
					BM_face_exists_overlap(vert_arr, 2, &f);
					convert_uv_to_new_face( v_buf1.orig_edge, f, &v1_u, &v1_v);
					convert_uv_to_new_face( v_buf2.orig_edge, f, &v2_u, &v2_v);
				}
			}
			face_index = BM_elem_index_get(f);
		}

		{
			float v1_uv[2] = { v1_u, v1_v };
			float v2_uv[2] = { v2_u, v2_v };

			printf("hej\n");
            bisect_search( v1_uv, v2_uv, eval, bm, e, face_index, cam_loc); 
		}

	}
}

static bool cusp_triangle(struct OpenSubdiv_EvaluatorDescr *eval, const float cam_loc[3], float co_arr[3][3],
		const bool b_arr[3], const float u_arr[3], const float v_arr[3], const int face_index, float cusp_co[3]){
	//If area is <= 1e-14, then we set the center of the triangle as the cusp
	//TODO the paper has 1e-20
	if( area_tri_v3(co_arr[0], co_arr[1], co_arr[2]) <= 1e-14 ){
		cent_tri_v3(cusp_co, co_arr[0], co_arr[1], co_arr[2]);
		return true;		
	}

	bool normal_sign_cross(const bool bool_arr[3]){
		int i;
		bool temp = bool_arr[0];
		for(i = 1; i < 3; i++){
			if(temp != bool_arr[i]){
				return true;
			}
		}
		return false;
	}

	bool k_r_sign_cross(const float u[3], const float v[3]){
		float k_r1 = get_k_r(eval, face_index, u[0], v[0], cam_loc);
		float k_r2 = get_k_r(eval, face_index, u[1], v[1], cam_loc);
		if( (k_r1 > 0) != (k_r2 > 0) ){
			//k_r sign crossing!
			return true;
		} else {
			//check last vert
			float k_r3 = get_k_r(eval, face_index, u[2], v[2], cam_loc);
			if( (k_r1 > 0) != (k_r3 > 0) ){
				return true;
			}
		}
		return false;
	}

    //print_m3("co_arr", co_arr);

	//printf("area: %.20f\n", area_tri_v3(co_arr[0], co_arr[1], co_arr[2]));

	{
		//Which edge is the longest in parameter space?
		float e1_len, e2_len, e3_len, du[3], dv[3];
		float new_co[3], uv_co[2];
		bool new_b;

		e1_len = len_v3v3(co_arr[0], co_arr[1]);
		e2_len = len_v3v3(co_arr[0], co_arr[2]);
		e3_len = len_v3v3(co_arr[1], co_arr[2]);

		/*
        printf("Before edge len\n");
		printf("e1: %f\n", e1_len);
		printf("e2: %f\n", e2_len);
		printf("e3: %f\n", e3_len);
        */

		//TODO there must be a way to reuse more code here...
		if(e1_len >= e2_len){
			if(e1_len >= e3_len){
				//e1
				float uv_1[] = {u_arr[0], v_arr[0]};
				float uv_2[] = {u_arr[1], v_arr[1]};
				
				interp_v2_v2v2( uv_co, uv_1, uv_2, 0.5f);
				openSubdiv_evaluateLimit(eval, face_index, uv_co[0], uv_co[1], new_co, du, dv);
				new_b = calc_if_B(cam_loc, new_co, du, dv);
				{
					bool new_b_arr[3] = { b_arr[0], new_b, b_arr[2] };
					float new_u_arr[3] = { u_arr[0], uv_co[0], u_arr[2] };
					float new_v_arr[3] = { v_arr[0], uv_co[1], v_arr[2] };

					if( normal_sign_cross(new_b_arr) && k_r_sign_cross(new_u_arr, new_v_arr) ){
						float new_co_arr[3][3];
						copy_m3_m3(new_co_arr, co_arr);
						copy_v3_v3(new_co_arr[1], new_co);
						if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp_co)){
							return true;
						}
					}
				}
				{
					bool new_b_arr[3] = { new_b, b_arr[1], b_arr[2] };
					float new_u_arr[3] = { uv_co[0], u_arr[1], u_arr[2] };
					float new_v_arr[3] = { uv_co[1], v_arr[1], v_arr[2] };

					if( normal_sign_cross(new_b_arr) && k_r_sign_cross(new_u_arr, new_v_arr) ){
						float new_co_arr[3][3];
						copy_m3_m3(new_co_arr, co_arr);
						copy_v3_v3(new_co_arr[0], new_co);
						if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp_co)){
							return true;
						}
					}
				}
			} else {
				//e3
				float uv_1[] = {u_arr[2], v_arr[2]};
				float uv_2[] = {u_arr[1], v_arr[1]};
				
				interp_v2_v2v2( uv_co, uv_1, uv_2, 0.5f);
				openSubdiv_evaluateLimit(eval, face_index, uv_co[0], uv_co[1], new_co, du, dv);
				new_b = calc_if_B(cam_loc, new_co, du, dv);

				{
					bool new_b_arr[3] = { b_arr[0], new_b, b_arr[2] };
					float new_u_arr[3] = { u_arr[0], uv_co[0], u_arr[2] };
					float new_v_arr[3] = { v_arr[0], uv_co[1], v_arr[2] };

					if( normal_sign_cross(new_b_arr) && k_r_sign_cross(new_u_arr, new_v_arr) ){
						float new_co_arr[3][3];
						copy_m3_m3(new_co_arr, co_arr);
						copy_v3_v3(new_co_arr[1], new_co);
						if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp_co)){
							return true;
						}
					}
				}
				{
					bool new_b_arr[3] = { b_arr[0], b_arr[1], new_b };
					float new_u_arr[3] = { u_arr[0], u_arr[1], uv_co[0] };
					float new_v_arr[3] = { v_arr[0], v_arr[1], uv_co[1] };

					if( normal_sign_cross(new_b_arr) && k_r_sign_cross(new_u_arr, new_v_arr) ){
						float new_co_arr[3][3];
						copy_m3_m3(new_co_arr, co_arr);
						copy_v3_v3(new_co_arr[2], new_co);
						if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp_co)){
							return true;
						}
					}
				}
			}
		} else if (e2_len >= e3_len) {
			//e2
			float uv_1[] = {u_arr[0], v_arr[0]};
			float uv_2[] = {u_arr[2], v_arr[2]};
				
			interp_v2_v2v2( uv_co, uv_1, uv_2, 0.5f);
			openSubdiv_evaluateLimit(eval, face_index, uv_co[0], uv_co[1], new_co, du, dv);
			new_b = calc_if_B(cam_loc, new_co, du, dv);

			{
				bool new_b_arr[3] = { new_b, b_arr[1], b_arr[2] };
				float new_u_arr[3] = { uv_co[0], u_arr[1], u_arr[2] };
				float new_v_arr[3] = { uv_co[1], v_arr[1], v_arr[2] };

				if( normal_sign_cross(new_b_arr) && k_r_sign_cross(new_u_arr, new_v_arr) ){
					float new_co_arr[3][3];
					copy_m3_m3(new_co_arr, co_arr);
					copy_v3_v3(new_co_arr[0], new_co);
					if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp_co)){
						return true;
					}
				}
			}
			{
				bool new_b_arr[3] = { b_arr[0], b_arr[1], new_b };
				float new_u_arr[3] = { u_arr[0], u_arr[1], uv_co[0] };
				float new_v_arr[3] = { v_arr[0], v_arr[1], uv_co[1] };

				if( normal_sign_cross(new_b_arr) && k_r_sign_cross(new_u_arr, new_v_arr) ){
					float new_co_arr[3][3];
					copy_m3_m3(new_co_arr, co_arr);
					copy_v3_v3(new_co_arr[2], new_co);
					if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp_co)){
						return true;
					}
				}
			}
		} else {
			//e3
			float uv_1[] = {u_arr[2], v_arr[2]};
			float uv_2[] = {u_arr[1], v_arr[1]};
				
			interp_v2_v2v2( uv_co, uv_1, uv_2, 0.5f);
			openSubdiv_evaluateLimit(eval, face_index, uv_co[0], uv_co[1], new_co, du, dv);
			new_b = calc_if_B(cam_loc, new_co, du, dv);

			{
				bool new_b_arr[3] = { b_arr[0], new_b, b_arr[2] };
				float new_u_arr[3] = { u_arr[0], uv_co[0], u_arr[2] };
				float new_v_arr[3] = { v_arr[0], uv_co[1], v_arr[2] };

				if( normal_sign_cross(new_b_arr) && k_r_sign_cross(new_u_arr, new_v_arr) ){
					float new_co_arr[3][3];
					copy_m3_m3(new_co_arr, co_arr);
					copy_v3_v3(new_co_arr[1], new_co);
					if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp_co)){
						return true;
					}
				}
			}
			{
				bool new_b_arr[3] = { b_arr[0], b_arr[1], new_b };
				float new_u_arr[3] = { u_arr[0], u_arr[1], uv_co[0] };
				float new_v_arr[3] = { v_arr[0], v_arr[1], uv_co[1] };

				if( normal_sign_cross(new_b_arr) && k_r_sign_cross(new_u_arr, new_v_arr) ){
					float new_co_arr[3][3];
					copy_m3_m3(new_co_arr, co_arr);
					copy_v3_v3(new_co_arr[2], new_co);
					if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp_co)){
						return true;
					}
				}
			}
		}

	}

	return false;
}

static void cusp_detection(BMesh *bm, BMesh *bm_orig, BLI_Buffer *new_vert_buffer,
						   BLI_Buffer *cusp_edges, struct OpenSubdiv_EvaluatorDescr *eval,
						   const float cam_loc[3]){
	BMFace *f;
	BMIter iter_f;

	int orig_verts = BM_mesh_elem_count(bm_orig, BM_VERT);
	int orig_edges = BM_mesh_elem_count(bm_orig, BM_EDGE);
    int initial_faces = BM_mesh_elem_count(bm, BM_FACE);
    int f_idx;

	BM_ITER_MESH_INDEX (f, &iter_f, bm, BM_FACES_OF_MESH, f_idx) {
		BMVert *vert;
		BMVert *vert_arr[3];
		BMIter iter_v;
		int vert_idx;
		bool first_vert, back_face, found_face, b_arr[3];
		first_vert = true;
		found_face = false;

        if( !(f_idx < initial_faces) ){
			//TODO perhaps insert every cusp edge after we iterated over all faces instead?
			//this might not be ok because will miss some faces that gets split
			break;
		}

		BM_ITER_ELEM_INDEX (vert, &iter_v, f, BM_VERTS_OF_FACE, vert_idx) {
			if(first_vert){
				first_vert = false;
				back_face = calc_if_B_nor(cam_loc, vert->co, vert->no);
				b_arr[vert_idx] = back_face;
			} else if (!found_face) {
				//If one or more of the verts do not have the same facing, then we want to look for cusps
				bool temp = calc_if_B_nor(cam_loc, vert->co, vert->no);
				b_arr[vert_idx] = temp;
				if(temp != back_face){
					found_face = true;
				}
			}
			vert_arr[vert_idx] = vert;
		}

		if(!found_face){
			continue;
		}

		//Find original mesh face + uv coords
		{
            float u_arr[3]; //array for u-coords (v1_u, v2_u ...)
			float v_arr[3];
			float co_arr[3][3];
			int i;

			BMEdge *edge_arr[] = {NULL, NULL, NULL};
			BMFace *orig_face = NULL;

            //check if all verts are on the orignal mesh
			for(i = 0; i < 3; i++){
				int v_idx = BM_elem_index_get(vert_arr[i]);
				
				//Copy coords for later use
				copy_v3_v3(co_arr[i], vert_arr[i]->co);

				if( (v_idx + 1) > orig_verts){
					Vert_buf v_buf = BLI_buffer_at(new_vert_buffer, Vert_buf, v_idx - orig_verts);
					u_arr[i] = v_buf.u;
					v_arr[i] = v_buf.v;
					if( v_buf.orig_face ){
						orig_face = v_buf.orig_face;
						vert_arr[i] = NULL;
					} else {
						//Make sure we don't pick one of the verts that we already have
						int idx1 = BM_elem_index_get(v_buf.orig_edge->v1); 
						
						if( (idx1 != BM_elem_index_get(vert_arr[ mod_i(i-1, 3) ]) ) && (idx1 != BM_elem_index_get(vert_arr[ mod_i(i+1, 3) ]) ) ){
							vert_arr[i] = v_buf.orig_edge->v1;							
						} else {
							//If v1 is a duplicate then v2 has to be unique
							vert_arr[i] = v_buf.orig_edge->v2;							
						}

						edge_arr[i] = v_buf.orig_edge;
					}
				} else {
					vert_arr[i] = BM_vert_at_index_find( bm_orig, v_idx );
				}
			}
            
			if(orig_face == NULL){
				BM_face_exists_overlap(vert_arr, 3, &orig_face);
			}

			for(i = 0; i < 3; i++){
				if(vert_arr[i] == NULL){
					continue;
				}

				if(edge_arr[i] != NULL){
					//Make use we have the correct uv coords
					convert_uv_to_new_face( edge_arr[i], orig_face, &u_arr[i], &v_arr[i]);
				} else {
					get_uv_coord(vert_arr[i], orig_face, &u_arr[i], &v_arr[i]);
				}
			}

			{
				int face_index = BM_elem_index_get(orig_face);
				//Check for k_r sign crossings
				float k_r1 = get_k_r(eval, face_index, u_arr[0], v_arr[0], cam_loc);
				float k_r2 = get_k_r(eval, face_index, u_arr[1], v_arr[1], cam_loc);
				bool k_r_crossing = false;
				if( (k_r1 > 0) != (k_r2 > 0) ){
					//k_r sign crossing!
					//printf("found k_r sign crossing\n");
					k_r_crossing = true;
				} else {
					//check last vert
					float k_r3 = get_k_r(eval, face_index, u_arr[2], v_arr[2], cam_loc);
					if( (k_r1 > 0) != (k_r3 > 0) ){
						//printf("found k_r sign crossing\n");
						k_r_crossing = true;
					}
				}

				if(k_r_crossing){
					float cusp_co[3];
					//Start looking for the cusp in the triangle
					if(cusp_triangle(eval, cam_loc, co_arr, b_arr, u_arr, v_arr, face_index, cusp_co)){
						//We found a cusp!
                        float a[3], b[3], c[3], d[3], no[3];
						float uv_1[2], uv_2[2];

						printf("Found a cusp point!\n");
						if(b_arr[0] == b_arr[1]){
							copy_v3_v3(a, co_arr[0]);
							copy_v3_v3(b, co_arr[1]);
							copy_v3_v3(c, co_arr[2]);

							uv_1[0] = u_arr[0];
							uv_2[0] = u_arr[1];

							uv_1[1] = v_arr[0];
							uv_2[1] = v_arr[1];
						} else if(b_arr[0] == b_arr[2]){
							copy_v3_v3(a, co_arr[0]);
							copy_v3_v3(b, co_arr[2]);
							copy_v3_v3(c, co_arr[1]);

							uv_1[0] = u_arr[0];
							uv_2[0] = u_arr[2];
                                   
							uv_1[1] = v_arr[0];
							uv_2[1] = v_arr[2];
						} else {
							copy_v3_v3(a, co_arr[1]);
							copy_v3_v3(b, co_arr[2]);
							copy_v3_v3(c, co_arr[0]);

							uv_1[0] = u_arr[1];
							uv_2[0] = u_arr[2];
                                   
							uv_1[1] = v_arr[1];
							uv_2[1] = v_arr[2];
						}
						sub_v3_v3v3(d, c, cusp_co);
						cross_v3_v3v3(no, f->no, d);
						
						if(isect_line_plane_v3(d, a, b, c, no)){
							float t, vec1[3], vec2[3];
							float edge_uv[2];
							BMEdge *edge;
							BMIter iter_e;
							sub_v3_v3v3(vec1, b, a);
							sub_v3_v3v3(vec2, d, a);

                            t = len_v3(vec2)/len_v3(vec1);

							BM_ITER_ELEM (edge, &iter_e, f, BM_EDGES_OF_FACE) {
								if( equals_v3v3(edge->v1->co, a) && equals_v3v3(edge->v2->co, b) ){
									//Found edge to split
									break;
								}
								if( equals_v3v3(edge->v1->co, b) && equals_v3v3(edge->v2->co, a) ){
									//Found edge to split
									t = fabsf(t - 1.0f);
									break;
								}
							}
							printf("t: %f\n", t);

							if(t > 1){
								print_m3("co_arr", co_arr);
								print_v3("cusp_co", cusp_co);
								//TODO fix this
								continue;
							}

							{
                                float P[3], du[3], dv[3];
								Vert_buf v_buf;
								int edge_idx = BM_elem_index_get(edge);
								int v1_idx = BM_elem_index_get(edge->v1);
								int v2_idx = BM_elem_index_get(edge->v2);

								interp_v2_v2v2(edge_uv, uv_1, uv_2, t);
								openSubdiv_evaluateLimit(eval, face_index, edge_uv[0], edge_uv[1], P, du, dv);

								if( len_v3v3(P, edge->v1->co) < 1e-3 || len_v3v3(P, edge->v2->co) < 1e-3 ){
									//Do not insert a new vert here
									continue;
								}

								if( (edge_idx < orig_edges) ){
									//Point on orig edge
									BMEdge *orig_e = BM_edge_at_index_find(bm_orig, edge_idx);
									v_buf.orig_edge = orig_e;
									v_buf.orig_face = NULL;
								} else if( edge_uv[0] == 0 || edge_uv[0] == 1 || edge_uv[1] == 0 || edge_uv[1] == 1 ){
									if( (v1_idx + 1) > orig_verts){
										Vert_buf v_buf_old = BLI_buffer_at(new_vert_buffer, Vert_buf, v1_idx - orig_verts);	
										v_buf.orig_edge = v_buf_old.orig_edge;
										v_buf.orig_face = NULL;
									} else if( (v2_idx + 1) > orig_verts) {
										Vert_buf v_buf_old = BLI_buffer_at(new_vert_buffer, Vert_buf, v2_idx - orig_verts);	
										v_buf.orig_edge = v_buf_old.orig_edge;
										v_buf.orig_face = NULL;
									}
								} else {
									//On orig face
									v_buf.orig_edge = NULL;
									v_buf.orig_face = orig_face;
								}

								v_buf.u = edge_uv[0];
								v_buf.v = edge_uv[1];

								BLI_buffer_append(new_vert_buffer, Vert_buf, v_buf);
								{
									int j;
									int cur_edge = BM_mesh_elem_count(bm, BM_EDGE);
                                    BMEdge *cusp_e = NULL;

									split_edge_and_move_vert(bm, edge, P, du, dv, edge_uv[0], edge_uv[1]);

                                    //Get the cusp edge
									for(j = 0; j < 3; j++){
										cur_edge++;
										cusp_e = BM_edge_at_index_find( bm, cur_edge );
										if( equals_v3v3(cusp_e->v1->co, c) && equals_v3v3(cusp_e->v2->co, P) ){
											break;
										}
										if( equals_v3v3(cusp_e->v1->co, P) && equals_v3v3(cusp_e->v2->co, c) ){
											break;
										}
										cusp_e = NULL;
									}
									if(cusp_e == NULL){
										printf("No cusp edge found!!!\n");
									} else {
										Cusp cusp;
										cusp.cusp_e = cusp_e;
										copy_v3_v3(cusp.cusp_co, cusp_co);
										BLI_buffer_append(cusp_edges, Cusp, cusp);
									}
								}
							}
						} else {
							//TODO handle this
							printf("ERROR: No intersection found!");
						}
					}
				}
			}
		}
	}                                                     
}

static void cusp_insertion(BMesh *bm, BLI_Buffer *cusp_edges){
	int cusp_i;

	for(cusp_i = 0; cusp_i < cusp_edges->count; cusp_i++){
		Cusp cusp = BLI_buffer_at(cusp_edges, Cusp, cusp_i);
		if( len_v3v3(cusp.cusp_co, cusp.cusp_e->v1->co) < 1e-2 || len_v3v3(cusp.cusp_co, cusp.cusp_e->v2->co) < 1e-2 ){
			//Do not insert a new vert here
			printf("skipped cups insert\n");
			continue;
		}
		split_edge_and_move_cusp(bm, cusp.cusp_e, cusp.cusp_co);
	}
}

static struct OpenSubdiv_EvaluatorDescr *create_osd_eval(BMesh *bm){
	//TODO create FAR meshes instead. (Perhaps the code in contour subdiv can help?)
	int subdiv_levels = 1;
	int no_of_verts = BM_mesh_elem_count(bm, BM_VERT);

	int j;
	int face_vert_index[4];
	float vert_array[3 * no_of_verts];

	BMIter iter_v, iter_f;
	BMVert *vert;
    BMFace *f;

	struct OpenSubdiv_EvaluatorDescr *osd_evaluator;
	osd_evaluator = openSubdiv_createEvaluatorDescr(no_of_verts);

	BM_ITER_MESH (f, &iter_f, bm, BM_FACES_OF_MESH) {
		BM_ITER_ELEM_INDEX (vert, &iter_v, f, BM_VERTS_OF_FACE, j) {
			face_vert_index[j] = BM_elem_index_get(vert);
		}
		//TODO this will only work with quad meshes. Set up checks so only quad meshes are passed through.
		openSubdiv_createEvaluatorDescrFace(osd_evaluator, 4, face_vert_index);
	}

	//TODO add check to see if this fails
	openSubdiv_finishEvaluatorDescr(osd_evaluator, subdiv_levels, OSD_SCHEME_CATMARK);

	BM_ITER_MESH_INDEX (vert, &iter_v, bm, BM_VERTS_OF_MESH, j) {
		vert_array[3*j] = vert->co[0];
		vert_array[3*j + 1] = vert->co[1];
		vert_array[3*j + 2] = vert->co[2];
	}

    openSubdiv_setEvaluatorCoarsePositions(osd_evaluator,
                                           vert_array,
                                           no_of_verts);

	return osd_evaluator;
}

static void debug_colorize(BMesh *bm, const float cam_loc[3]){
	BMIter iter;
	BMFace *f;
	float P[3];

	BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH){
		BM_face_calc_center_mean(f, P);
		if( calc_if_B_nor(cam_loc, P, f->no) ){
			f->mat_nr = 1;
		}
	}
}

/* bmesh only function */
static DerivedMesh *mybmesh_do(DerivedMesh *dm, MyBMeshModifierData *mmd, float cam_loc[3])
{

	DerivedMesh *result;
	BMesh *bm_orig, *bm;

	//CCGSubSurf *ss;
    struct OpenSubdiv_EvaluatorDescr *osd_eval;

	bm = DM_to_bmesh(dm, true);
	//Keep a copy of the original mesh
	bm_orig = DM_to_bmesh(dm, true);
    
    osd_eval = create_osd_eval(bm);

	// (6.1) Initialization
    verts_to_limit(bm, osd_eval);

	BM_mesh_triangulate(bm, MOD_TRIANGULATE_QUAD_FIXED, MOD_TRIANGULATE_NGON_BEAUTY, false, NULL, NULL);

	if( mmd->camera_ob == NULL){
		//Can't proceed without camera obj
		result = CDDM_from_bmesh(bm, true);
		BM_mesh_free(bm);
		BM_mesh_free(bm_orig);
		openSubdiv_deleteEvaluatorDescr(osd_eval);
		return result;
	}
	{
		BLI_buffer_declare_static(Vert_buf, new_vert_buffer, BLI_BUFFER_NOP, 32);
		BLI_buffer_declare_static(Cusp, cusp_edges, BLI_BUFFER_NOP, 32);

		split_BB_FF_edges(bm, bm_orig, osd_eval, cam_loc, &new_vert_buffer);

		// (6.2) Contour Insertion

        //TODO implement vertex shift (as an alternative to edge split)

        cusp_detection(bm, bm_orig, &new_vert_buffer, &cusp_edges, osd_eval, cam_loc);

		contour_insertion(bm, bm_orig, &new_vert_buffer, &cusp_edges, osd_eval, cam_loc);

        cusp_insertion(bm, &cusp_edges);

		debug_colorize(bm, cam_loc);
		BLI_buffer_free(&new_vert_buffer);
		BLI_buffer_free(&cusp_edges);
	}
	result = CDDM_from_bmesh(bm, true);

	BM_mesh_free(bm);
	BM_mesh_free(bm_orig);

	openSubdiv_deleteEvaluatorDescr(osd_eval);

	result->dirty |= DM_DIRTY_NORMALS;

	//TODO use this to check if we need to subdivide the mesh to get a quad mesh.
	if(0){
	BMIter iter;
	BMFace *f, *f_next;
	int sides = 4;
	/* use the mutable iterator so we can remove data as its looped over */
	BM_ITER_MESH_MUTABLE (f, f_next, &iter, bm, BM_FACES_OF_MESH) {
		if (f->len == sides) {
			BM_face_kill(bm, f);
		}
	}
	}

	return result;
}

/* MyBMesh */
static void initData(ModifierData *md)
{
	MyBMeshModifierData *mmd = (MyBMeshModifierData *) md;

	mmd->camera_ob = NULL;
}

static void copyData(ModifierData *md, ModifierData *target)
{
	/* NOTE: you might want to copy some settings manually (if simple memcpy isn't enough)*/
#if 0
	MyBMeshModifierData *mmd  = (MyBMeshModifierData *)md;
	MyBMeshModifierData *tmmd = (MyBMeshModifierData *)target;
#endif
	modifier_copyData_generic(md, target);
}

static DerivedMesh *applyModifier(ModifierData *md, Object *ob,
                                  DerivedMesh *dm,
                                  ModifierApplyFlag UNUSED(flag))
{
	DerivedMesh *result;
	float cam_loc[3];
	MyBMeshModifierData *mmd = (MyBMeshModifierData *)md;

	if(mmd->camera_ob){

		copy_v3_v3(cam_loc, mmd->camera_ob->loc);
		printf("Cam loc:\n");
		printf("1: %f\n", cam_loc[0]);
		printf("2: %f\n", cam_loc[1]);
		printf("3: %f\n", cam_loc[2]);
		//convert camera origin from world coord to the modifier obj local coords
		mul_m4_v3(ob->obmat, cam_loc);
		printf("Cam loc 2:\n");
		printf("1: %f\n", cam_loc[0]);
		printf("2: %f\n", cam_loc[1]);
		printf("3: %f\n", cam_loc[2]);
	}

	if (!(result = mybmesh_do(dm, mmd, cam_loc))) {
		return dm;
	}

	return result;
}

static bool dependsOnNormals(ModifierData *UNUSED(md))
{
	return false;
}

ModifierTypeInfo modifierType_MyBMesh = {
	/* name */              "MyBMesh",
	/* structName */        "MyBMeshModifierData",
	/* structSize */        sizeof(MyBMeshModifierData),
    /* type */              eModifierTypeType_Constructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	                        eModifierTypeFlag_SupportsEditmode |
	                        eModifierTypeFlag_EnableInEditmode,

	/* copyData */          copyData,
	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     applyModifier,
	/* applyModifierEM */   NULL,
	/* initData */          initData,
	/* requiredDataMask */  NULL,
	/* freeData */          NULL,
	/* isDisabled */        NULL,
	//TODO implement depgraph
	/* updateDepgraph */    NULL,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */  dependsOnNormals,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};
