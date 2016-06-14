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

#include "bmesh.h"
#include "bmesh_tools.h"

#include "MOD_util.h"

#include "depsgraph_private.h"

//TODO this modifier depends on OSD. So if it's not compiled in, remove this modifier
#include <opensubdiv/osdutil/evaluator_capi.h>

struct OpenSubdiv_EvaluatorDescr;

typedef struct {
	BMVert *vert;
	BMEdge *orig_edge;
	BMFace *orig_face;
	float u, v;
} Vert_buf;

typedef struct {
	BMEdge *cusp_e;
	BMFace *orig_face;
	float cusp_co[3];
	float cusp_no[3];

	float u, v;
} Cusp;

typedef struct {
	BMFace *face;
	//Should be front or back facing?
	bool back_f;
} IncoFace;

typedef struct {
	BMesh *bm;
	BMesh *bm_orig;

	float cam_loc[3];

	BLI_Buffer *new_vert_buffer;
	BLI_Buffer *shifted_verts;
	BLI_Buffer *cusp_edges;
	BLI_Buffer *C_verts;
	//Radial edge vert start idx
	int radi_start_idx;

	//are we in the cusp insertion step
	bool is_cusp;

	struct OpenSubdiv_EvaluatorDescr *eval;
} MeshData;

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
					case 1 :
						u = 1, v = 0;
						break;
					case 2 :
						u = v = 1;
						break;
					case 3 :
						u = 0, v = 1;
						break;
					default:
						u = v = 0;
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

static BMVert *split_edge_and_move_vert(BMesh *bm, BMEdge *edge, const float new_pos[3],
									const float du[3], const float dv[3]){
	//Split edge one time and move the created vert to new_pos

	BMVert *vert, *temp_v;
	BMOperator div_op;
	BMOIter oiter;

	printf("Split edge!\n");

	BM_mesh_elem_hflag_disable_all(bm, BM_EDGE, BM_ELEM_TAG, false);
	BM_elem_flag_enable(edge, BM_ELEM_TAG);
	BMO_op_initf(bm, &div_op, BMO_FLAG_DEFAULTS,
			"subdivide_edges edges=%he cuts=%i quad_corner_type=%i use_single_edge=%b",
			BM_ELEM_TAG, 1, SUBD_STRAIGHT_CUT, true, "geom.out");

	BMO_op_exec(bm, &div_op);

	//Get the newly created vertex
	BMO_ITER (temp_v, &oiter, div_op.slots_out, "geom.out", BM_VERT) {
		vert = temp_v;
	}

	copy_v3_v3(vert->co, new_pos);
	//Adjust vert normal to the limit normal
	cross_v3_v3v3(vert->no, dv, du);
	normalize_v3(vert->no);

	BMO_op_finish(bm, &div_op);
	return vert;
}


static BMVert* split_edge_and_move_cusp(BMesh *bm, BMEdge *edge, const float new_pos[3], const float new_no[3]){
	//Split edge one time and move the created vert to new_pos

	BMVert *vert, *temp_v;
	BMOperator div_op;
	BMOIter oiter;
	printf("Cusp split!\n");

	BM_mesh_elem_hflag_disable_all(bm, BM_EDGE, BM_ELEM_TAG, false);
	BM_elem_flag_enable(edge, BM_ELEM_TAG);
	BMO_op_initf(bm, &div_op, BMO_FLAG_DEFAULTS,
			"subdivide_edges edges=%he cuts=%i quad_corner_type=%i use_single_edge=%b",
			BM_ELEM_TAG, 1, SUBD_STRAIGHT_CUT, true, "geom.out");

	BMO_op_exec(bm, &div_op);

	//Get the newly created vertex
	BMO_ITER (temp_v, &oiter, div_op.slots_out, "geom.out", BM_VERT) {
		vert = temp_v;
	}

	copy_v3_v3(vert->co, new_pos);
	copy_v3_v3(vert->no, new_no);

	BMO_op_finish(bm, &div_op);
	return vert;
}

static bool get_uv_coord(BMVert *vert, BMFace *f, float *u, float *v){
	//Get U,V coords of a vertex
	int i;
	BMIter iter;
	BMVert *vert_iter;

	BM_ITER_ELEM_INDEX (vert_iter, &iter, f, BM_VERTS_OF_FACE, i) {
		if(vert == vert_iter){
			switch(i){
				case 1 :
					*u = 1, *v = 0;
					break;
				case 2 :
					*u = *v = 1;
					break;
				case 3 :
					*u = 0, *v = 1;
					break;
				default:
					*u = *v = 0;
					break;
			}
			return true;
		}
	}

	return false;
}

static void split_BB_FF_edges(MeshData *m_d) {
	//Split BB,FF edges if they have sign crossings

	int i, face_index;
	BMIter iter_f, iter_e;
	BMEdge *e;
	BMFace *f;
	BMVert *v1, *v2;
	float v1_u, v1_v, v2_u, v2_v;
	bool is_B;
	int orig_edges = BM_mesh_elem_count(m_d->bm_orig, BM_EDGE);
	int initial_edges = BM_mesh_elem_count(m_d->bm, BM_EDGE);

	//Do 10 samples but don't check end and start point
	float step = 1.0f/11.0f;
	float step_arr[] = { step*5.0f, step*6.0f, step*4.0f, step*7.0f, step*3.0f,
		step*8.0f, step*2.0f, step*9.0f, step*1.0f, step*10.0f };

	BM_ITER_MESH_INDEX (e, &iter_e, m_d->bm, BM_EDGES_OF_MESH, i) {
		Vert_buf v_buf;

		if( !(i < initial_edges) ){
			//Now we are working on edges we added in this function
			break;
		}

		//TODO perhaps we need to use the limit surface normal

		is_B = calc_if_B_nor(m_d->cam_loc, e->v1->co, e->v1->no);

		if( is_B  != calc_if_B_nor(m_d->cam_loc, e->v2->co, e->v2->no) ){
			//This is not a FF or BB edge
			continue;
		}

		if( i < orig_edges ){
			//This edge exists on the original mesh
			//TODO why do I have to use find? Segfault otherwise...
			//remember to replace the rest of "at_index"
			//why is table_ensure not fixing the assert?
			BMEdge *orig_e = BM_edge_at_index_find(m_d->bm_orig, i);

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
			v_buf.orig_face = f;
		} else {
			BMVert *vert_arr[2];
			bool found_face;

			//This should be safe because the vert count is still the same as the original mesh.
			v1 = BM_vert_at_index_find( m_d->bm_orig, BM_elem_index_get( e->v1 ) );
			v2 = BM_vert_at_index_find( m_d->bm_orig, BM_elem_index_get( e->v2 ) );

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
			openSubdiv_evaluateLimit(m_d->eval, face_index, v1_u, v1_v, P1, du, dv);

			is_B = calc_if_B(m_d->cam_loc, P1, du, dv);

			openSubdiv_evaluateLimit(m_d->eval, face_index, v2_u, v2_v, P2, du, dv);

			if( is_B  != calc_if_B(m_d->cam_loc, P2, du, dv) ){
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

				for(i=0; i < 10; i++){
					v = step_arr[i];
					openSubdiv_evaluateLimit(m_d->eval, face_index, u, v, P, du, dv);
					if( calc_if_B(m_d->cam_loc, P, du, dv) != is_B ){
						split_edge_and_move_vert(m_d->bm, e, P, du, dv);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
						break;
					}
				}
			} else if ( v1_v == v2_v ){
				v = v1_v;

				for(i=0; i < 10; i++){
					u = step_arr[i];
					openSubdiv_evaluateLimit(m_d->eval, face_index, u, v, P, du, dv);
					if( calc_if_B(m_d->cam_loc, P, du, dv) != is_B ){
						split_edge_and_move_vert(m_d->bm, e, P, du, dv);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
						break;
					}
				}
			} else {
				bool alt_diag;
				if((v1_u == 0 && v1_v == 0) || (v2_u == 0 && v2_v == 0)){
					alt_diag = false;
				} else {
					alt_diag = true;
				}
				for(i=0; i < 10; i++){
					if(alt_diag){
						u = 1.0f - step_arr[i];
					} else {
						u = step_arr[i];
					}

					v = step_arr[i];
					openSubdiv_evaluateLimit(m_d->eval, face_index, u, v, P, du, dv);
					if( calc_if_B(m_d->cam_loc, P, du, dv) != is_B ){
						split_edge_and_move_vert(m_d->bm, e, P, du, dv);
						v_buf.u = u;
						v_buf.v = v;
						BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
						break;
					}
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

static void convert_uv_to_new_face(BMEdge *e, BMFace *old_f, BMFace *f, float *u, float *v){
	//convert the old u/v coords to the new face coord

	float v1_u, v1_v, v2_u, v2_v;
	float old_v1_u, old_v1_v, old_v2_u, old_v2_v;

	if(old_f == f){
		//Already have the correct uv coords
		return;
	}

	get_uv_coord(e->v1, f, &v1_u, &v1_v);
	get_uv_coord(e->v2, f, &v2_u, &v2_v);

	get_uv_coord(e->v1, old_f, &old_v1_u, &old_v1_v);
	get_uv_coord(e->v2, old_f, &old_v2_u, &old_v2_v);

	//Which axis are we moving along?
	if( *u == 0.0f || *u == 1.0f ){
		//Along v axis
		if( v1_u == v2_u ){
			//Still along the v axis in the new face
			if(v1_v != old_v1_v){
				*v = 1.0f - *v;
			}
			*u = v1_u;
		} else {
			//Changed axis to u
			if(v1_u != old_v1_v){
				*u = 1.0f - *v;
			} else {
				*u = *v;
			}
			*v = v1_v;
		}
	} else {
		//Along u axis
		if( v1_v == v2_v ){
			//Still along the u axis in the new face
			if(v1_u != old_v1_u){
				*u = 1.0f - *u;
			}
			*v = v1_v;
		} else {
			//Changed axis to v
			if(v1_v != old_v1_u){
				*v = 1.0f - *u;
			} else {
				*v = *u;
			}
			*u = v1_u;
		}
	}

}

static bool append_vert(BLI_Buffer *C_verts, BMVert *vert){
	int vert_i;

	//check if vert is already in the buffer
	for(vert_i = 0; vert_i < C_verts->count; vert_i++){
		if( vert == BLI_buffer_at(C_verts, BMVert*, vert_i)){
			return false;
		}
	}
	BLI_buffer_append(C_verts, BMVert*, vert);
	return true;
}

static Vert_buf* get_shift_vert( BMVert *vert, MeshData *m_d ){
	int vert_i;

	//check if vert is in the buffer
	for(vert_i = 0; vert_i < m_d->shifted_verts->count; vert_i++){
		Vert_buf *buf = &BLI_buffer_at(m_d->shifted_verts, Vert_buf, vert_i);
		if( vert == buf->vert ){
			return buf;
		}
	}
	return NULL;
}

static void add_shifted_vert( BMVert *vert, BMFace *orig_face, float uv[2], MeshData *m_d ){
	Vert_buf *buf = get_shift_vert( vert, m_d );

	if(buf != NULL){
		buf->orig_face = orig_face;
		buf->u = uv[0];
		buf->v = uv[1];
	} else {
		Vert_buf new_buf;
		new_buf.orig_face = orig_face;
		new_buf.vert = vert;
		new_buf.u = uv[0];
		new_buf.v = uv[1];
		BLI_buffer_append(m_d->shifted_verts, Vert_buf, new_buf);
	}
}

static bool check_and_shift(BMVert *vert, const float new_loc[3], const float new_no[3], MeshData *m_d){
	//TODO add all shiftability checks from the paper
	typedef struct {
		float no[3];
	} Normal;

	float old_loc[3];

	copy_v3_v3( old_loc, vert->co );

	// Will the shift create folds?
	// TODO perhaps only checking for huge normal changes in not enough?
	{

		BLI_buffer_declare_static(Normal, old_normals, BLI_BUFFER_NOP, 32);
		BMFace* f;
		BMIter iter_f;
		int i = 0;

		//Copy old face normals
		BM_ITER_ELEM (f, &iter_f, vert, BM_FACES_OF_VERT) {
			Normal nor;
			BM_face_calc_normal(f, nor.no);
			BLI_buffer_append(&old_normals, Normal, nor);
		}

		copy_v3_v3(vert->co, new_loc);

		//Check if the new position changed any of the normals drastically (potential fold)
		BM_ITER_ELEM (f, &iter_f, vert, BM_FACES_OF_VERT) {
			float no[3];
			float old_no[3];

			BM_face_calc_normal(f, no);

			copy_v3_v3( old_no, BLI_buffer_at(&old_normals, Normal, i).no );
			if( dot_v3v3( old_no, no ) < 0.5f ){
				//Big change in normal dir, potential fold, abort
				copy_v3_v3(vert->co, old_loc);
				printf("Skipped shift vert!\n");
				BLI_buffer_free(&old_normals);
				return false;
			}

			i++;
		}

		//Move the vert back for future checks
		copy_v3_v3(vert->co, old_loc);

		BLI_buffer_free(&old_normals);
	}

	//Will this shift a cusp edge
	if( !m_d->is_cusp ){
		int edge_i;
		BMEdge *edge;
		BMIter iter_e;

		BM_ITER_ELEM (edge, &iter_e, vert, BM_EDGES_OF_VERT) {
			for(edge_i = 0; edge_i < m_d->cusp_edges->count; edge_i++){
				BMEdge* cusp_edge = BLI_buffer_at(m_d->cusp_edges, Cusp, edge_i).cusp_e;
				if( edge == cusp_edge ){
					return false;
				}
			}
		}
	}

	//Check if the shift might/will cause a CCC face
	{
		BLI_buffer_declare_static(BMVert*, c_verts, BLI_BUFFER_NOP, 32);
		BMFace* f;
		BMIter iter_f;

		BM_ITER_ELEM (f, &iter_f, vert, BM_FACES_OF_VERT) {
			BMEdge* e;
			BMIter iter_e;
			int num_cross = 0; //number of zero crossing edges
			BM_ITER_ELEM (e, &iter_e, f, BM_EDGES_OF_FACE) {
				if( e->v1 != vert && e->v2 != vert ){
					int vert_i;
					int edge_c_verts = 0;
					for(vert_i = 0; vert_i < m_d->C_verts->count; vert_i++){
						BMVert* C_vert = BLI_buffer_at(m_d->C_verts, BMVert*, vert_i);
						if( e->v1 == C_vert ) {
							edge_c_verts++;

							if( append_vert( &c_verts, e->v1 ) ){
								num_cross++;
							}
						}
						if( e->v2 == C_vert ){
							edge_c_verts++;

							if( append_vert( &c_verts, e->v2 ) ){
								num_cross++;
							}
						}
					}

					if ( edge_c_verts >= 2 ){
						//TODO if > 2 then we have added duplicate verts to C_verts, should not happen!
						BLI_buffer_free(&c_verts);
						return false;
					}

					if( edge_c_verts == 0 && calc_if_B_nor(m_d->cam_loc, e->v1->co, e->v1->no) != calc_if_B_nor(m_d->cam_loc, e->v2->co, e->v2->no) ){
						num_cross++;
					}

					if( num_cross > 2 ) {
						BLI_buffer_free(&c_verts);
						return false;
					}

				}
			}
		}
		BLI_buffer_free(&c_verts);
	}
	//Adjust vert normal to the limit normal
	copy_v3_v3(vert->no, new_no);
	copy_v3_v3(vert->co, new_loc);
	return true;
}

static void mult_face_search( BMFace *f, BMFace *f2, const float v1_uv[2], const float v2_uv[2], BMEdge *e, MeshData *m_d ){
	//Try to find a vert that is connected to both faces
	BMVert *vert;
	BMFace *face;
	BMIter iter_f, iter_v;
	bool found_vert = false;

	BM_ITER_ELEM (vert, &iter_v, f, BM_VERTS_OF_FACE) {
		if( !BM_vert_is_boundary(vert) && BM_vert_edge_count(vert) == 4 && BM_vert_face_count(vert) == 4 ){
			BM_ITER_ELEM (face, &iter_f, vert, BM_FACES_OF_VERT) {
				if(face == f2){
					found_vert = true;
					break;
				}
			}
		}
		if(found_vert){
			break;
		}
	}

	if( !found_vert ){
		//We can't easily interpolate this edge, do not try to insert a new vertex here
		printf("Couldn't find any suitable interpolation vertex!\n");
		return;
	}

	//The vert we found is the center of our new 2d coordinate system "ab".
	//Convert the uv coords to ab coords
	{
		//UV axis is the four quadrats, "+" is the selected vertex
		//   b (up)
		//   ^
		// |2|1|
		// --+--> a (right)
		// |3|4|

		// The uv points in the origin, up, right corner of the face
		float face_uv[4][3][2];
		float ab_start[2], ab_end[2];
		int edge_idx = 0, face_idx = 0;
		BMEdge *up, *down, *left, *right, *cur_edge;
		BMFace *quadrants[4];

		BMLoop *first_loop = BM_face_vert_share_loop( vert->e->l->f, vert );
		BMLoop *cur_loop = first_loop;

		cur_edge = vert->e;

		do {
			switch(edge_idx){
				case 0 :
					right = cur_edge;
					break;
				case 1 :
					up = cur_edge;
					break;
				case 2 :
					left = cur_edge;
					break;
				default:
					down = cur_edge;
					break;
			}
			edge_idx++;
		} while (((cur_loop = BM_vert_step_fan_loop(cur_loop, &cur_edge)) != first_loop) && (cur_loop != NULL));

		do {
			float u, v;
			face = cur_loop->f;

			//Save face for later use
			quadrants[face_idx] = face;

			get_uv_coord(vert, face, &u, &v);
			face_uv[ face_idx ][0][0] = u;
			face_uv[ face_idx ][0][1] = v;


			get_uv_coord(cur_loop->next->v, face, &u, &v);

			// The BMLoop edge consist of the current vertex and the vertex in cur_loop->next->v
			if( cur_loop->e == right ){
				face_uv[ face_idx ][2][0] = u;
				face_uv[ face_idx ][2][1] = v;
			} else if( cur_loop->e == up ){
				face_uv[ face_idx ][1][0] = u;
				face_uv[ face_idx ][1][1] = v;
			} else if( cur_loop->e == left ){
				face_uv[ face_idx ][2][0] = u;
				face_uv[ face_idx ][2][1] = v;
			} else {
				//down
				face_uv[ face_idx ][1][0] = u;
				face_uv[ face_idx ][1][1] = v;
			}

			get_uv_coord(cur_loop->prev->v, face, &u, &v);
			if( cur_loop->prev->e == right ){
				face_uv[ face_idx ][2][0] = u;
				face_uv[ face_idx ][2][1] = v;
			} else if( cur_loop->prev->e == up ){
				face_uv[ face_idx ][1][0] = u;
				face_uv[ face_idx ][1][1] = v;
			} else if( cur_loop->prev->e == left ){
				face_uv[ face_idx ][2][0] = u;
				face_uv[ face_idx ][2][1] = v;
			} else {
				//down
				face_uv[ face_idx ][1][0] = u;
				face_uv[ face_idx ][1][1] = v;
			}

			//Convert the supplied uv coords to ab coords
			if( face == f || face == f2 ){
				float a, b;

				if( face == f ){
					u = v1_uv[0];
					v = v1_uv[1];
				} else {
					u = v2_uv[0];
					v = v2_uv[1];
				}

				//Convert the u axis
				if( face_uv[ face_idx ][0][0] == 0 ){
					if( face_uv[ face_idx ][2][0] == 0 ){
						//b coord is mapped to u
						b = u;
					} else {
						//a coord is mapped to u
						a = u;
					}
				} else {
					if( face_uv[ face_idx ][2][0] == 0 ){
						//a coord is mapped to u
						a = 1.0f - u;
					} else {
						//b coord is mapped to u
						b = 1.0f - u;
					}
				}

				//Convert the v axis
				if( face_uv[ face_idx ][0][1] == 0 ){
					if( face_uv[ face_idx ][1][1] == 0 ){
						//a coord is mapped to v
						a = v;
					} else {
						//b coord is mapped to v
						b = v;
					}
				} else {
					if( face_uv[ face_idx ][1][1] == 0 ){
						//b coord is mapped to v
						b = 1.0f - v;
					} else {
						//a coord is mapped to v
						a = 1.0f - v;
					}
				}

				switch( face_idx ){
					case 1 :
						//2nd quadrant
						a = -a;
						break;
					case 2 :
						//3nd quadrant
						a = -a;
						b = -b;
						break;
					case 3 :
						//4th quadrant
						b = -b;
						break;
					default :
						//1st quadrant
						break;
				}

				//first or second face?
				if( face == f ){
					ab_start[0] = a;
					ab_start[1] = b;
				} else {
					ab_end[0] = a;
					ab_end[1] = b;
				}

			}
			face_idx++;
		} while (((cur_loop = BM_vert_step_fan_loop(cur_loop, &cur_edge)) != first_loop) && (cur_loop != NULL));

		//Now we can begin interpolating along the edge
		{
			float face_dir, uv_P[2], uv_1[2], uv_2[2],  P[3], du[3], dv[3], new_no[3];
			float step = 0.5f;
			float step_len = 0.25f;
			float cur_ab[2];
			int i, face_index, q_idx;
			float v1_face = get_facing_dir_nor(m_d->cam_loc, e->v1->co, e->v1->no);
			BMFace *cur_face;

			for( i = 0; i < 10; i++ ){
				interp_v2_v2v2( cur_ab, ab_start, ab_end, step);
				if( cur_ab[0] < 0 ){
					if( cur_ab[1] < 0 ){
						q_idx = 2;
					} else {
						q_idx = 1;
					}
				} else {
					if( cur_ab[1] < 0 ){
						q_idx = 3;
					} else {
						q_idx = 0;
					}
				}

				cur_face = quadrants[q_idx];

				interp_v2_v2v2( uv_1, face_uv[q_idx][0], face_uv[q_idx][2], fabs(cur_ab[0]) );
				interp_v2_v2v2( uv_2, face_uv[q_idx][0], face_uv[q_idx][1], fabs(cur_ab[1]) );

				if( face_uv[q_idx][0][0] == 1 ){
					if( face_uv[q_idx][2][0] == 1 ){
						uv_1[0] = 0;
					}
					if( face_uv[q_idx][1][0] == 1 ){
						uv_2[0] = 0;
					}
				}

				if( face_uv[q_idx][0][1] == 1 ){
					if( face_uv[q_idx][2][1] == 1 ){
						uv_1[1] = 0;
					}
					if( face_uv[q_idx][1][1] == 1 ){
						uv_2[1] = 0;
					}
				}

				add_v2_v2v2( uv_P, uv_1, uv_2 );

				face_index = BM_elem_index_get(cur_face);
				openSubdiv_evaluateLimit(m_d->eval, face_index, uv_P[0], uv_P[1], P, du, dv);

				face_dir = get_facing_dir(m_d->cam_loc, P, du, dv);

				if( fabs(face_dir) < 1e-20 ){
					//We got lucky and found the zero crossing!
					printf("--->> got lucky\n");
					break;
				}

				if( (face_dir < 0) == (v1_face < 0) ){
					step += step_len;
				} else {
					step -= step_len;
				}
				step_len = step_len/2.0f;
			}

			cross_v3_v3v3(new_no, dv, du);
			normalize_v3(new_no);

			if( len_v3v3(P, e->v1->co) < BM_edge_calc_length(e) * 0.2f ){
				if(check_and_shift(e->v1, P, new_no, m_d) ){
					//Do not insert a new vert here, shift it instead
					append_vert(m_d->C_verts, e->v1);
					add_shifted_vert( e->v1, cur_face, uv_P, m_d );
					return;
				} else if( len_v3v3(P, e->v1->co) < BM_edge_calc_length(e) * 0.01f ){
					append_vert(m_d->C_verts, e->v1);
					return;
				}
			} else if (len_v3v3(P, e->v2->co) < BM_edge_calc_length(e) * 0.2f ){
				if(check_and_shift(e->v2, P, new_no, m_d) ){
				//Do not insert a new vert here, shift it instead
				append_vert(m_d->C_verts, e->v2);
				add_shifted_vert( e->v2, cur_face, uv_P, m_d );
				return;
				} else if( len_v3v3(P, e->v2->co) < BM_edge_calc_length(e) * 0.01f ){
					append_vert(m_d->C_verts, e->v2);
					return;
				}
			}

			{
				//Insert a new vert

				Vert_buf new_buf;

				BMVert *vert = split_edge_and_move_vert(m_d->bm, e, P, du, dv);
				append_vert(m_d->C_verts, vert);

				new_buf.orig_face = cur_face;
				new_buf.orig_edge = NULL;
				new_buf.u = uv_P[0];
				new_buf.v = uv_P[1];
				BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, new_buf);
			}
		}
	}

}

static bool bisect_search(const float v1_uv[2], const float v2_uv[2], BMEdge *e, BMFace *orig_face, float uv_result[2], MeshData *m_d ){
	//Search edge for sign crossing and split it!
	int i;
	float face_dir, uv_P[2], P[3], du[3], dv[3], new_no[3];
	float step = 0.5f;
	float step_len = 0.25f;
	float v1_face = get_facing_dir_nor(m_d->cam_loc, e->v1->co, e->v1->no);
	int face_index = BM_elem_index_get(orig_face);

	for( i = 0; i < 10; i++){
		interp_v2_v2v2( uv_P, v1_uv, v2_uv, step);
		openSubdiv_evaluateLimit(m_d->eval, face_index, uv_P[0], uv_P[1], P, du, dv);
		face_dir = get_facing_dir(m_d->cam_loc, P, du, dv);

		if( fabs(face_dir) < 1e-20 ){
			//We got lucky and found the zero crossing!
			printf("--->> got lucky\n");
			break;
		}

		if( (face_dir < 0) == (v1_face < 0) ){
			step += step_len;
		} else {
			step -= step_len;
		}
		step_len = step_len/2.0f;
	}

	cross_v3_v3v3(new_no, dv, du);
	normalize_v3(new_no);

	if( len_v3v3(P, e->v1->co) < BM_edge_calc_length(e) * 0.2f ){
		if( check_and_shift(e->v1, P, new_no, m_d) ){
			//Do not insert a new vert here, shift it instead
			append_vert(m_d->C_verts, e->v1);
			add_shifted_vert( e->v1, orig_face, uv_P, m_d );
			return false;
		} else if( len_v3v3(P, e->v1->co) < BM_edge_calc_length(e) * 0.01f ){
			append_vert(m_d->C_verts, e->v1);
			return false;
		}
	} else if (len_v3v3(P, e->v2->co) < BM_edge_calc_length(e) * 0.2f ){
		if( check_and_shift(e->v2, P, new_no, m_d) ){
			//Do not insert a new vert here, shift it instead
			append_vert(m_d->C_verts, e->v2);
			add_shifted_vert( e->v2, orig_face, uv_P, m_d );
			return false;
		} else if( len_v3v3(P, e->v2->co) < BM_edge_calc_length(e) * 0.01f ){
			append_vert(m_d->C_verts, e->v2);
			return false;
		}
	}

	copy_v2_v2(uv_result, uv_P);
	{
		BMVert *vert = split_edge_and_move_vert(m_d->bm, e, P, du, dv);
		append_vert(m_d->C_verts, vert);
	}
	return true;
}

static void search_edge( const int i, BMEdge *e, MeshData *m_d){

	float v1_u, v1_v, v2_u, v2_v;
	int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
	int orig_edges = BM_mesh_elem_count(m_d->bm_orig, BM_EDGE);
	int v1_idx = BM_elem_index_get(e->v1);
	int v2_idx = BM_elem_index_get(e->v2);
	Vert_buf v_buf1, v_buf2;
	BMFace *f, *f2;
	BMEdge *orig_e = NULL;
	BMVert *v1 = NULL, *v2 = NULL;
	bool v1_has_face = false, v2_has_face = false, diff_faces = false;

	if( (v1_idx + 1) > orig_verts){
		v_buf1 = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v1_idx - orig_verts);
		v1_u = v_buf1.u;
		v1_v = v_buf1.v;
		if( v_buf1.orig_edge == NULL ){
			v1_has_face = true;
		}
	} else {
		v1 = BM_vert_at_index_find( m_d->bm_orig, v1_idx );
	}
	if( (v2_idx + 1) > orig_verts){
		v_buf2 = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v2_idx - orig_verts);
		v2_u = v_buf2.u;
		v2_v = v_buf2.v;
		if( v_buf2.orig_edge == NULL ){
			v2_has_face = true;
		}
	} else {
		v2 = BM_vert_at_index_find( m_d->bm_orig, v2_idx );
	}

	if( v1 && v2 ){

		if( i < orig_edges ){
			//this edge is on the original mesh
			BMIter iter_f;

			//TODO why do I have to use find? Segfault otherwise...
			//remember to replace the rest of "at_index"
			//why is table_ensure not fixing the assert?
			orig_e = BM_edge_at_index_find(m_d->bm_orig, i);

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
		if( v2_has_face ){
			f = v_buf2.orig_face;
		} else {
			BMVert *vert_arr[3];

			vert_arr[0] = v_buf2.orig_edge->v1;
			vert_arr[1] = v_buf2.orig_edge->v2;
			vert_arr[2] = v1;
			//TODO check if get face fails
			if( v_buf2.orig_edge->v1 == v2 || v_buf2.orig_edge->v2 == v2 ){
				BM_face_exists_overlap(vert_arr, 2, &f);
			} else {
				BM_face_exists_overlap(vert_arr, 3, &f);
			}
			convert_uv_to_new_face( v_buf2.orig_edge, v_buf2.orig_face, f, &v2_u, &v2_v);
		}
		get_uv_coord(v1, f, &v1_u, &v1_v);
	} else if ( v2 ){
		if( v1_has_face ){
			f = v_buf1.orig_face;
		} else {
			BMVert *vert_arr[3];

			vert_arr[0] = v_buf1.orig_edge->v1;
			vert_arr[1] = v_buf1.orig_edge->v2;
			vert_arr[2] = v2;
			//TODO check if get face fails
			if( v_buf1.orig_edge->v1 == v2 || v_buf1.orig_edge->v2 == v2 ){
				BM_face_exists_overlap(vert_arr, 2, &f);
			} else {
				BM_face_exists_overlap(vert_arr, 3, &f);
			}
			convert_uv_to_new_face( v_buf1.orig_edge, v_buf1.orig_face, f, &v1_u, &v1_v);
		}
		get_uv_coord(v2, f, &v2_u, &v2_v);
	} else {
		if( v1_has_face || v2_has_face ){
			if( v1_has_face && v2_has_face ){
				if( v_buf1.orig_face != v_buf2.orig_face ){
					diff_faces = true;
					f2 = v_buf2.orig_face;
				}
				f = v_buf1.orig_face;
			} else if ( v1_has_face ) {
				f = v_buf1.orig_face;
				convert_uv_to_new_face( v_buf2.orig_edge, v_buf2.orig_face, f, &v2_u, &v2_v);
			} else {
				f = v_buf2.orig_face;
				convert_uv_to_new_face( v_buf1.orig_edge, v_buf1.orig_face, f, &v1_u, &v2_v);
			}
		} else {
			//No orig face. So this in on a orig edge. So just get the face from the v1 edge
			BMVert *vert_arr[] = {v_buf1.orig_edge->v1 ,v_buf1.orig_edge->v2};
			BM_face_exists_overlap(vert_arr, 2, &f);
			convert_uv_to_new_face( v_buf1.orig_edge, v_buf1.orig_face, f, &v1_u, &v1_v);
			convert_uv_to_new_face( v_buf2.orig_edge, v_buf2.orig_face, f, &v2_u, &v2_v);
		}
	}

	{
		//TODO rewrite the above checks so that we do not need to do this check
		//Check if one of the verts has been shifted or not
		Vert_buf* vert1 = get_shift_vert( e->v1, m_d );
		Vert_buf* vert2 = get_shift_vert( e->v2, m_d );

		if( vert1 != NULL ){
			if(vert1->orig_face != f){
				if(diff_faces){
					f = vert1->orig_face;
				} else {
					diff_faces = true;
					f2 = f;
					f = vert1->orig_face;
				}
			}
			v1_u = vert1->u;
			v1_v = vert1->v;
		}

		if( vert2 != NULL ){
			if(diff_faces){
				if(vert2->orig_face != f2){
					f2 = vert2->orig_face;
				}
			} else if(vert2->orig_face != f) {
				diff_faces = true;
				f2 = vert2->orig_face;
			}

			v2_u = vert2->u;
			v2_v = vert2->v;
		}

		//Check if a shifted vert caused the verts to be on the same face
		if(diff_faces && f == f2){
			diff_faces = false;
		}
	}

	{
		Vert_buf new_buf;
		float uv_result[2];
		float v1_uv[2] = { v1_u, v1_v };
		float v2_uv[2] = { v2_u, v2_v };

		if(diff_faces) {
			//The edge spawns over multiple original edges, try to interpolate along this edge.
			//If it fails, do not insert any new verts here
			printf("Mult face search\n");
			mult_face_search( f, f2, v1_uv, v2_uv, e, m_d );
			return;
		}

		if( (v1_u == 0 && v2_u == 0) || (v1_u == 1 && v2_u == 1) ||
			(v1_v == 0 && v2_v == 0) || (v1_v == 1 && v2_v == 1) )
		{
			//Along an original edge, save orig face for uv conversion
			new_buf.orig_face = f;
			if( v1 && v2 ){
				new_buf.orig_edge = orig_e;
			} else if ( v1 ){
				new_buf.orig_edge = v_buf2.orig_edge;
			} else {
				new_buf.orig_edge = v_buf1.orig_edge;
			}

		} else {
			new_buf.orig_face = f;
			new_buf.orig_edge = NULL;
		}

		if( bisect_search( v1_uv, v2_uv, e, f, uv_result, m_d) ){
			//if a new vert is inserted add it to the buffer
			new_buf.u = uv_result[0];
			new_buf.v = uv_result[1];
			BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, new_buf);
		}
	}
}

static void contour_insertion( MeshData *m_d ) {
	int i, cusp_i;
	BMEdge *e;
	BMIter iter_e;

	int initial_edges = BM_mesh_elem_count(m_d->bm, BM_EDGE);

	printf("Buffer count: %d\n", m_d->new_vert_buffer->count);

	BM_ITER_MESH_INDEX (e, &iter_e, m_d->bm, BM_EDGES_OF_MESH, i) {
		if( !(i < initial_edges) ){
			//Now we are working on edges we added in this function
			break;
		}

		{
			bool cont = false;
			//Check if this is a cusp edge
			for(cusp_i = 0; cusp_i < m_d->cusp_edges->count; cusp_i++){
				Cusp cusp = BLI_buffer_at(m_d->cusp_edges, Cusp, cusp_i);
				if(cusp.cusp_e == e){
					//Do not split this edge yet
					printf("skipped cusp edge\n");
					cont = true;
					break;
				}
			}
			if( cont ){
				continue;
			}
		}

		//TODO perhaps we need to use the limit surface normal
		if( calc_if_B_nor(m_d->cam_loc, e->v1->co, e->v1->no) == calc_if_B_nor(m_d->cam_loc, e->v2->co, e->v2->no) ){
			//This is not a FB or BF edge
			continue;
		}

		//Check if the edge already has a C vert
		{
			int vert_i;
			bool skip_edge = false;
			for(vert_i = 0; vert_i < m_d->C_verts->count; vert_i++){
				BMVert* C_vert = BLI_buffer_at(m_d->C_verts, BMVert*, vert_i);
				if( e->v1 == C_vert || e->v2 == C_vert){
					skip_edge = true;
				}
			}

			if(skip_edge) {
				continue;
			}
		}
		search_edge( i, e, m_d );
	}
}

static bool cusp_triangle(struct OpenSubdiv_EvaluatorDescr *eval, const float cam_loc[3], float co_arr[3][3],
		const bool b_arr[3], const float u_arr[3], const float v_arr[3], const int face_index, Cusp *cusp){
	//If area is <= 1e-14, then we set the center of the triangle as the cusp
	//TODO the paper has 1e-20
	if( area_tri_v3(co_arr[0], co_arr[1], co_arr[2]) <= 1e-14 ){
		float cusp_co[3];
		float cusp_no[3];
		float du[3], dv[3];

		float uv1[3] = { u_arr[0], v_arr[0], 0 };
		float uv2[3] = { u_arr[1], v_arr[1], 0 };
		float uv3[3] = { u_arr[2], v_arr[2], 0 };

		float res_uv[3];

		cent_tri_v3(res_uv, uv1, uv2, uv3);
		openSubdiv_evaluateLimit(eval, face_index, res_uv[0], res_uv[1], cusp_co, du, dv);
		cross_v3_v3v3(cusp_no, dv, du);

		copy_v3_v3(cusp->cusp_co, cusp_co);
		copy_v3_v3(cusp->cusp_no, cusp_no);
		cusp->u = res_uv[0];
		cusp->v = res_uv[1];
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
						if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp)){
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
						if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp)){
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
						if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp)){
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
						if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp)){
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
					if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp)){
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
					if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp)){
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
					if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp)){
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
					if(cusp_triangle(eval, cam_loc, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, cusp)){
						return true;
					}
				}
			}
		}

	}

	return false;
}

static BMFace *get_orig_face(int orig_verts, const BMVert *vert_arr_in[3], float u_arr[3], float v_arr[3], float co_arr[3][3], MeshData *m_d){
	int i;

	BMEdge *edge_arr[] = {NULL, NULL, NULL};
	BMEdge *edge_face_arr[] = {NULL, NULL, NULL};
	BMFace *orig_face = NULL;
	BMVert *vert_arr[3] = {vert_arr_in[0], vert_arr_in[1], vert_arr_in[2]};

	//check if all verts are on the orignal mesh
	for(i = 0; i < 3; i++){
		int v_idx = BM_elem_index_get(vert_arr[i]);

		//Copy coords for later use
		copy_v3_v3(co_arr[i], vert_arr[i]->co);

		if( (v_idx + 1) > orig_verts){
			Vert_buf v_buf = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v_idx - orig_verts);
			u_arr[i] = v_buf.u;
			v_arr[i] = v_buf.v;
			if( v_buf.orig_edge == NULL ){
				orig_face = v_buf.orig_face;
				vert_arr[i] = NULL;
			} else {
				//Make sure we don't pick one of the verts that we already have
				int idx1 = BM_elem_index_get(v_buf.orig_edge->v1);
				BMVert *temp_v1 = vert_arr[ mod_i(i-1, 3) ];
				BMVert *temp_v2 = vert_arr[ mod_i(i+1, 3) ];

				if( ( temp_v1 != NULL && idx1 != BM_elem_index_get(temp_v1) ) &&
						( temp_v2 != NULL && idx1 != BM_elem_index_get(temp_v2) ) )
				{
					vert_arr[i] = v_buf.orig_edge->v1;
				} else {
					//If v1 is a duplicate then v2 has to be unique
					vert_arr[i] = v_buf.orig_edge->v2;
				}

				edge_arr[i] = v_buf.orig_edge;
				edge_face_arr[i] = v_buf.orig_face;
			}
		} else {
			vert_arr[i] = BM_vert_at_index_find( m_d->bm_orig, v_idx );
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
			convert_uv_to_new_face( edge_arr[i], edge_face_arr[i], orig_face, &u_arr[i], &v_arr[i]);
		} else {
			get_uv_coord(vert_arr[i], orig_face, &u_arr[i], &v_arr[i]);
		}
	}
	return orig_face;
}

static void cusp_detection( MeshData *m_d ){
	BMFace *f;
	BMIter iter_f;

	int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
	int orig_edges = BM_mesh_elem_count(m_d->bm_orig, BM_EDGE);
	int initial_faces = BM_mesh_elem_count(m_d->bm, BM_FACE);
	int f_idx;

	BM_ITER_MESH_INDEX (f, &iter_f, m_d->bm, BM_FACES_OF_MESH, f_idx) {
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
				back_face = calc_if_B_nor(m_d->cam_loc, vert->co, vert->no);
				b_arr[vert_idx] = back_face;
			} else if (!found_face) {
				//If one or more of the verts do not have the same facing, then we want to look for cusps
				bool temp = calc_if_B_nor(m_d->cam_loc, vert->co, vert->no);
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
			BMFace *orig_face = get_orig_face(orig_verts, vert_arr, u_arr, v_arr, co_arr, m_d);

			{
				int face_index = BM_elem_index_get(orig_face);
				//Check for k_r sign crossings
				float k_r1 = get_k_r(m_d->eval, face_index, u_arr[0], v_arr[0], m_d->cam_loc);
				float k_r2 = get_k_r(m_d->eval, face_index, u_arr[1], v_arr[1], m_d->cam_loc);
				bool k_r_crossing = false;
				if( (k_r1 > 0) != (k_r2 > 0) ){
					//k_r sign crossing!
					//printf("found k_r sign crossing\n");
					k_r_crossing = true;
				} else {
					//check last vert
					float k_r3 = get_k_r(m_d->eval, face_index, u_arr[2], v_arr[2], m_d->cam_loc);
					if( (k_r1 > 0) != (k_r3 > 0) ){
						//printf("found k_r sign crossing\n");
						k_r_crossing = true;
					}
				}

				if(k_r_crossing){
					Cusp cusp;
					cusp.orig_face = orig_face;
					//Start looking for the cusp in the triangle
					if(cusp_triangle(m_d->eval, m_d->cam_loc, co_arr, b_arr, u_arr, v_arr, face_index, &cusp)){
						//We found a cusp!
						float uv_1[2], uv_2[2], uv_3[2];
						BMEdge *edge;
						BMVert *cusp_e_vert;

						printf("Found a cusp point!\n");
						if(b_arr[0] == b_arr[1]){
							uv_1[0] = u_arr[0];
							uv_2[0] = u_arr[1];
							uv_3[0] = u_arr[2];

							uv_1[1] = v_arr[0];
							uv_2[1] = v_arr[1];
							uv_3[1] = v_arr[2];
							edge = BM_edge_exists( vert_arr[0], vert_arr[1] );
							cusp_e_vert = vert_arr[2];
						} else if(b_arr[0] == b_arr[2]){
							uv_1[0] = u_arr[0];
							uv_2[0] = u_arr[2];
							uv_3[0] = u_arr[1];

							uv_1[1] = v_arr[0];
							uv_2[1] = v_arr[2];
							uv_3[1] = v_arr[1];
							edge = BM_edge_exists( vert_arr[0], vert_arr[2] );
							cusp_e_vert = vert_arr[1];
						} else {
							uv_1[0] = u_arr[1];
							uv_2[0] = u_arr[2];
							uv_3[0] = u_arr[0];

							uv_1[1] = v_arr[1];
							uv_2[1] = v_arr[2];
							uv_3[1] = v_arr[0];
							edge = BM_edge_exists( vert_arr[1], vert_arr[2] );
							cusp_e_vert = vert_arr[0];
						}

						{
							float P[3], du[3], dv[3];
							float edge_uv[2];
							float cusp_uv[2] = {cusp.u, cusp.v};
							Vert_buf v_buf;
							int edge_idx = BM_elem_index_get(edge);
							int v1_idx = BM_elem_index_get(edge->v1);
							int v2_idx = BM_elem_index_get(edge->v2);

							if( isect_line_line_v2_point( uv_1, uv_2, uv_3, cusp_uv, edge_uv ) != ISECT_LINE_LINE_CROSS ){
								printf("Couldn't find intersection point to edge from cusp!\n");
								//TODO this is a big error so quit instead
								continue;
							}

							openSubdiv_evaluateLimit(m_d->eval, face_index, edge_uv[0], edge_uv[1], P, du, dv);

							//Check if we should use an existing edge (no new verts)
							if( len_v3v3(P, edge->v1->co) < BM_edge_calc_length(edge) * 0.2f ||
									len_v3v3(P, edge->v2->co) < BM_edge_calc_length(edge) * 0.2f ){

								BMVert *edge_vert;
								BMEdge *cusp_edge;
								BMIter iter_e;
								float new_no[3];

								if( len_v3v3(P, edge->v1->co) < BM_edge_calc_length(edge) * 0.2f ){
									edge_vert = edge->v1;
								} else {
									edge_vert = edge->v2;
								}

								BM_ITER_ELEM (cusp_edge, &iter_e, f, BM_EDGES_OF_FACE) {
									if( cusp_edge != edge && (cusp_edge->v1 == edge_vert || cusp_edge->v2 == edge_vert) ){
										//Found edge
										break;
									}
								}

								cross_v3_v3v3(new_no, dv, du);
								normalize_v3(new_no);

								//Can we shift this vertex?
								if( check_and_shift(edge_vert, P, new_no, m_d) ){
									cusp.cusp_e = cusp_edge;
									add_shifted_vert( edge_vert , orig_face, edge_uv, m_d );
									BLI_buffer_append(m_d->cusp_edges, Cusp, cusp);

									printf("Used existing edge for cusp!\n");
									continue;
								}
							}

							if( (edge_idx < orig_edges) ){
								//Point on orig edge
								BMEdge *orig_e = BM_edge_at_index_find(m_d->bm_orig, edge_idx);
								v_buf.orig_edge = orig_e;
								v_buf.orig_face = orig_face;
							} else if( edge_uv[0] == 0 || edge_uv[0] == 1 || edge_uv[1] == 0 || edge_uv[1] == 1 ){
								if( (v1_idx + 1) > orig_verts){
									Vert_buf v_buf_old = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v1_idx - orig_verts);
									v_buf.orig_edge = v_buf_old.orig_edge;
									v_buf.orig_face = orig_face;
								} else if( (v2_idx + 1) > orig_verts) {
									Vert_buf v_buf_old = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v2_idx - orig_verts);
									v_buf.orig_edge = v_buf_old.orig_edge;
									v_buf.orig_face = orig_face;
								}
							} else {
								//On orig face
								v_buf.orig_edge = NULL;
								v_buf.orig_face = orig_face;
							}

							v_buf.u = edge_uv[0];
							v_buf.v = edge_uv[1];

							BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
							{
								int j;
								int cur_edge = BM_mesh_elem_count(m_d->bm, BM_EDGE);
								BMEdge *cusp_e = NULL;

								split_edge_and_move_vert(m_d->bm, edge, P, du, dv);

								//Get the cusp edge
								for(j = 0; j < 3; j++){
									cur_edge++;
									cusp_e = BM_edge_at_index_find( m_d->bm, cur_edge );
									if( cusp_e->v1 == cusp_e_vert || cusp_e->v2 == cusp_e_vert ){
										break;
									}
									cusp_e = NULL;
								}
								if(cusp_e == NULL){
									printf("No cusp edge found!!!\n");
								} else {
									cusp.cusp_e = cusp_e;
									BLI_buffer_append(m_d->cusp_edges, Cusp, cusp);
								}
							}
						}
					}
				}
			}
		}
	}
}

static void cusp_insertion(MeshData *m_d){
	int cusp_i;

	m_d->is_cusp = true;

	for(cusp_i = 0; cusp_i < m_d->cusp_edges->count; cusp_i++){
		BMVert *vert;
		Vert_buf new_buf;
		Cusp cusp = BLI_buffer_at(m_d->cusp_edges, Cusp, cusp_i);

		if( len_v3v3(cusp.cusp_co, cusp.cusp_e->v1->co) < BM_edge_calc_length(cusp.cusp_e) * 0.2f ){
			if(  BM_vert_edge_count(cusp.cusp_e->v1) == 4 && check_and_shift(cusp.cusp_e->v1, cusp.cusp_co, cusp.cusp_no, m_d) ){
				float uv_P[2] = { cusp.u, cusp.v };
				append_vert(m_d->C_verts, cusp.cusp_e->v1);
				add_shifted_vert( cusp.cusp_e->v1, cusp.orig_face, uv_P, m_d );
				continue;
			}
		} else if( len_v3v3(cusp.cusp_co, cusp.cusp_e->v2->co) < BM_edge_calc_length(cusp.cusp_e) * 0.2f ){
			if(  BM_vert_edge_count(cusp.cusp_e->v2) == 4 && check_and_shift(cusp.cusp_e->v2, cusp.cusp_co, cusp.cusp_no, m_d) ){
				float uv_P[2] = { cusp.u, cusp.v };
				append_vert(m_d->C_verts, cusp.cusp_e->v2);
				add_shifted_vert( cusp.cusp_e->v2, cusp.orig_face, uv_P, m_d );
				continue;
			}
		}

		vert = split_edge_and_move_cusp(m_d->bm, cusp.cusp_e, cusp.cusp_co, cusp.cusp_no);
		append_vert(m_d->C_verts, vert);

		new_buf.orig_face = cusp.orig_face;
		new_buf.orig_edge = NULL;
		new_buf.u = cusp.u;
		new_buf.v = cusp.v;

		BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, new_buf);
	}
	m_d->is_cusp = false;
}

static bool rad_triangle(struct OpenSubdiv_EvaluatorDescr *eval, const float rad_plane_no1[3], const float rad_plane_no2[3],
		const float CC1_pos[3], const float CC2_pos[3], float co_arr[3][3],
		const bool b_arr[3], const float u_arr[3], const float v_arr[3], const int face_index, float uv[2]){
	//If area is <= 1e-14, then we set the center of the triangle as the cusp
	//TODO the paper has 1e-20
	if( area_tri_v3(co_arr[0], co_arr[1], co_arr[2]) <= 1e-14 ){
		uv[0] = ( u_arr[0] + u_arr[1] + u_arr[3] ) / 3.0f;
		uv[1] = ( v_arr[0] + v_arr[1] + v_arr[3] ) / 3.0f;

		return true;
	}

	bool rad1_sign_cross(const bool bool_arr[3]){
		int i;
		bool temp = bool_arr[0];
		for(i = 1; i < 3; i++){
			if(temp != bool_arr[i]){
				return true;
			}
		}
		return false;
	}

	bool rad2_sign_cross(const float u[3], const float v[3]){
		float temp[3], P[3], val1, val2;

		openSubdiv_evaluateLimit(eval, face_index, u[0], v[0], P, NULL, NULL);

		sub_v3_v3v3(temp, P, CC2_pos);
		val1 = dot_v3v3(rad_plane_no2, temp);

		openSubdiv_evaluateLimit(eval, face_index, u[1], v[1], P, NULL, NULL);

		sub_v3_v3v3(temp, P, CC2_pos);
		val2 = dot_v3v3(rad_plane_no2, temp);

		if( (val1 > 0) != (val2 > 0) ){
			//rad sign crossing!
			return true;
		} else {
			//check last vert
			float val3;

			openSubdiv_evaluateLimit(eval, face_index, u[2], v[2], P, NULL, NULL);

			sub_v3_v3v3(temp, P, CC2_pos);
			val3 = dot_v3v3(rad_plane_no2, temp);
			if( (val1 > 0) != (val3 > 0) ){
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
		float new_co[3], temp[3], uv_co[2];
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
				sub_v3_v3v3(temp, new_co, CC1_pos);
				new_b = (dot_v3v3(rad_plane_no1, temp) > 0);
				{
					bool new_b_arr[3] = { b_arr[0], new_b, b_arr[2] };
					float new_u_arr[3] = { u_arr[0], uv_co[0], u_arr[2] };
					float new_v_arr[3] = { v_arr[0], uv_co[1], v_arr[2] };

					if( rad1_sign_cross(new_b_arr) && rad2_sign_cross(new_u_arr, new_v_arr) ){
						float new_co_arr[3][3];
						copy_m3_m3(new_co_arr, co_arr);
						copy_v3_v3(new_co_arr[1], new_co);
						if(rad_triangle(eval, rad_plane_no1, rad_plane_no2, CC1_pos, CC2_pos, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, uv)){
							return true;
						}
					}
				}
				{
					bool new_b_arr[3] = { new_b, b_arr[1], b_arr[2] };
					float new_u_arr[3] = { uv_co[0], u_arr[1], u_arr[2] };
					float new_v_arr[3] = { uv_co[1], v_arr[1], v_arr[2] };

					if( rad1_sign_cross(new_b_arr) && rad2_sign_cross(new_u_arr, new_v_arr) ){
						float new_co_arr[3][3];
						copy_m3_m3(new_co_arr, co_arr);
						copy_v3_v3(new_co_arr[0], new_co);
						if(rad_triangle(eval, rad_plane_no1, rad_plane_no2, CC1_pos, CC2_pos, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, uv)){
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
				sub_v3_v3v3(temp, new_co, CC1_pos);
				new_b = (dot_v3v3(rad_plane_no1, temp) > 0);

				{
					bool new_b_arr[3] = { b_arr[0], new_b, b_arr[2] };
					float new_u_arr[3] = { u_arr[0], uv_co[0], u_arr[2] };
					float new_v_arr[3] = { v_arr[0], uv_co[1], v_arr[2] };

					if( rad1_sign_cross(new_b_arr) && rad2_sign_cross(new_u_arr, new_v_arr) ){
						float new_co_arr[3][3];
						copy_m3_m3(new_co_arr, co_arr);
						copy_v3_v3(new_co_arr[1], new_co);
						if(rad_triangle(eval, rad_plane_no1, rad_plane_no2, CC1_pos, CC2_pos, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, uv)){
							return true;
						}
					}
				}
				{
					bool new_b_arr[3] = { b_arr[0], b_arr[1], new_b };
					float new_u_arr[3] = { u_arr[0], u_arr[1], uv_co[0] };
					float new_v_arr[3] = { v_arr[0], v_arr[1], uv_co[1] };

					if( rad1_sign_cross(new_b_arr) && rad2_sign_cross(new_u_arr, new_v_arr) ){
						float new_co_arr[3][3];
						copy_m3_m3(new_co_arr, co_arr);
						copy_v3_v3(new_co_arr[2], new_co);
						if(rad_triangle(eval, rad_plane_no1, rad_plane_no2, CC1_pos, CC2_pos, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, uv)){
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
			sub_v3_v3v3(temp, new_co, CC1_pos);
			new_b = (dot_v3v3(rad_plane_no1, temp) > 0);

			{
				bool new_b_arr[3] = { new_b, b_arr[1], b_arr[2] };
				float new_u_arr[3] = { uv_co[0], u_arr[1], u_arr[2] };
				float new_v_arr[3] = { uv_co[1], v_arr[1], v_arr[2] };

				if( rad1_sign_cross(new_b_arr) && rad2_sign_cross(new_u_arr, new_v_arr) ){
					float new_co_arr[3][3];
					copy_m3_m3(new_co_arr, co_arr);
					copy_v3_v3(new_co_arr[0], new_co);
					if(rad_triangle(eval, rad_plane_no1, rad_plane_no2, CC1_pos, CC2_pos, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, uv)){
						return true;
					}
				}
			}
			{
				bool new_b_arr[3] = { b_arr[0], b_arr[1], new_b };
				float new_u_arr[3] = { u_arr[0], u_arr[1], uv_co[0] };
				float new_v_arr[3] = { v_arr[0], v_arr[1], uv_co[1] };

				if( rad1_sign_cross(new_b_arr) && rad2_sign_cross(new_u_arr, new_v_arr) ){
					float new_co_arr[3][3];
					copy_m3_m3(new_co_arr, co_arr);
					copy_v3_v3(new_co_arr[2], new_co);
					if(rad_triangle(eval, rad_plane_no1, rad_plane_no2, CC1_pos, CC2_pos, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, uv)){
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
			sub_v3_v3v3(temp, new_co, CC1_pos);
			new_b = (dot_v3v3(rad_plane_no1, temp) > 0);

			{
				bool new_b_arr[3] = { b_arr[0], new_b, b_arr[2] };
				float new_u_arr[3] = { u_arr[0], uv_co[0], u_arr[2] };
				float new_v_arr[3] = { v_arr[0], uv_co[1], v_arr[2] };

				if( rad1_sign_cross(new_b_arr) && rad2_sign_cross(new_u_arr, new_v_arr) ){
					float new_co_arr[3][3];
					copy_m3_m3(new_co_arr, co_arr);
					copy_v3_v3(new_co_arr[1], new_co);
					if(rad_triangle(eval, rad_plane_no1, rad_plane_no2, CC1_pos, CC2_pos, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, uv)){
						return true;
					}
				}
			}
			{
				bool new_b_arr[3] = { b_arr[0], b_arr[1], new_b };
				float new_u_arr[3] = { u_arr[0], u_arr[1], uv_co[0] };
				float new_v_arr[3] = { v_arr[0], v_arr[1], uv_co[1] };

				if( rad1_sign_cross(new_b_arr) && rad2_sign_cross(new_u_arr, new_v_arr) ){
					float new_co_arr[3][3];
					copy_m3_m3(new_co_arr, co_arr);
					copy_v3_v3(new_co_arr[2], new_co);
					if(rad_triangle(eval, rad_plane_no1, rad_plane_no2, CC1_pos, CC2_pos, new_co_arr, new_b_arr, new_u_arr, new_v_arr, face_index, uv)){
						return true;
					}
				}
			}
		}

	}

	return false;
}


static bool is_C_vert(BMVert *v, BLI_Buffer *C_verts){
	int vert_j;
	for(vert_j = 0; vert_j < C_verts->count; vert_j++){
		BMVert *c_vert = BLI_buffer_at(C_verts, BMVert*, vert_j);
		if( c_vert == v ){
			return true;
		}
	}
	return false;
}

static bool point_inside_v2(const float mat[3][3], const float point[2], BMFace *f){
	//TODO maybe add a sanity check to see if the face is not a quad or a triangle
	float mat_coords[f->len][2];
	BMVert *vert;
	BMIter iter_v;
	int vert_idx;

	BM_ITER_ELEM_INDEX (vert, &iter_v, f, BM_VERTS_OF_FACE, vert_idx) {
		mul_v2_m3v3(mat_coords[vert_idx], mat, vert->co);
	}

	if( f->len == 3 ){
		return isect_point_tri_v2(point, mat_coords[0], mat_coords[1], mat_coords[2]);
	}
	return isect_point_quad_v2(point, mat_coords[0], mat_coords[1], mat_coords[2], mat_coords[3]);
}

static bool point_inside(const float mat[3][3], const float point[3], BMFace *f){
	float mat_new_pos[2];
	mul_v2_m3v3(mat_new_pos, mat, point);

	return point_inside_v2(mat, mat_new_pos, f);
}

static void get_uv_point(BMFace *face, float uv[2], const float point_v2[2], const float mat[3][3] ){
	int vert_idx;
	float st[4][2];

	BMVert *v;
	BMIter iter_v;

	BM_ITER_ELEM_INDEX (v, &iter_v, face, BM_VERTS_OF_FACE, vert_idx) {
		switch(vert_idx){
			case 1 :
				mul_v2_m3v3(st[1], mat, v->co);
				break;
			case 2 :
				mul_v2_m3v3(st[2], mat, v->co);
				break;
			case 3 :
				mul_v2_m3v3(st[3], mat, v->co);
				break;
			default:
				mul_v2_m3v3(st[0], mat, v->co);
				break;
		}
	}

	resolve_quad_uv_v2(uv, point_v2, st[0], st[1], st[2], st[3]);
}

static bool poke_and_move(BMFace *f, const float new_pos[3], const float du[3], const float dv[3], MeshData *m_d){
	BMVert *vert, *temp_v;
	BMOperator poke_op;
	BMOIter oiter;

	BMEdge *edge = NULL;
	bool rot_edge = false;
	float mat[3][3];
	float new_norm[3];

	cross_v3_v3v3(new_norm, dv, du);
	normalize_v3(new_norm);

	axis_dominant_v3_to_m3(mat, new_norm);

	// BM_face_point_inside_test is too inaccurate to use here as some overhangs are missed with it.
	// TODO check if point_inside_test needs to be replaced in other places also
	if( !point_inside(mat, new_pos, f) ){
		BMIter iter_e;
		BMIter iter_f;
		BMEdge *e;
		BMFace *face;

		rot_edge = true;

		BM_ITER_ELEM (e, &iter_e, f, BM_EDGES_OF_FACE){
			BM_ITER_ELEM (face, &iter_f, e, BM_FACES_OF_EDGE){
				if( face != f ){

					if( point_inside(mat, new_pos, face) ){
						edge = e;
						break;
					}
				}
			}
		}
	}

	if( rot_edge ){
		if( edge == NULL || !BM_edge_rotate_check(edge) ){
			//Do not insert a radial edge here
			return false;
		}
	}

	BM_mesh_elem_hflag_disable_all(m_d->bm, BM_FACE, BM_ELEM_TAG, false);
	BM_elem_flag_enable(f, BM_ELEM_TAG);

	BMO_op_initf(m_d->bm, &poke_op, BMO_FLAG_DEFAULTS,
			"poke faces=%hf", BM_ELEM_TAG, "verts.out");

	BMO_op_exec(m_d->bm, &poke_op);

	//Get the newly created vertex
	BMO_ITER (temp_v, &oiter, poke_op.slots_out, "verts.out", BM_VERT) {
		vert = temp_v;
	}

	copy_v3_v3(vert->co, new_pos);

	//Adjust vert normal to the limit normal
	copy_v3_v3(vert->no, new_norm);

	if( rot_edge ){
		BM_edge_rotate(m_d->bm, edge, true, 0);
		printf("rotated edge!\n");
	}


	BMO_op_finish(m_d->bm, &poke_op);

	return true;
}

static void mult_radi_search( BMFace *diff_f[3], const float cent[3], const float edge1_mid[3], const float edge2_mid[3],
							const float val_1, const float val_2,
							const float rad_plane_no[3], const float C_vert_pos[3], BMFace *poke_face, MeshData *m_d ){
	//Try to find a vert that is connected to both faces
	BMVert *vert;
	BMFace *face;
	BMIter iter_f, iter_v;
	int edge_count, f_idx;
	bool found_vert = false;
	float mat[3][3];

	for(f_idx = 0; f_idx < 3; f_idx++){
		BMFace *f = diff_f[f_idx];
		BM_ITER_ELEM (vert, &iter_v, f, BM_VERTS_OF_FACE) {
			if( !BM_vert_is_boundary(vert) && BM_vert_edge_count(vert) == BM_vert_face_count(vert) ){
				bool e1 = false;
				bool e2 = false;
				bool f_cent = false;

				axis_dominant_v3_to_m3(mat, vert->no);

				BM_ITER_ELEM (face, &iter_f, vert, BM_FACES_OF_VERT) {
					if( point_inside(mat, edge1_mid, face) ){
						e1 = true;
					}
					if( point_inside(mat, edge2_mid, face) ){
						e2 = true;
					}
					if( point_inside(mat, cent, face) ){
						f_cent = true;
					}
					if(e1 && e2 && f_cent){
						edge_count = BM_vert_edge_count(vert);
						found_vert = true;
						break;
					}
				}
			}
			if(found_vert){
				break;
			}
		}
		if(found_vert){
			break;
		}
	}

	if( !found_vert ){
		//We can't easily interpolate this edge, do not try to insert a new vertex here
		printf("Couldn't find any suitable interpolation vertex!\n");
		return;
	}

	{
		int edge_idx = 0;

		BMLoop *first_loop = BM_face_vert_share_loop( vert->e->l->f, vert );
		BMLoop *cur_loop = first_loop;
		BMEdge *cur_edge;
		BMFace *faces[edge_count];

		cur_edge = vert->e;

		do {
			faces[edge_idx] = cur_loop->f;
			edge_idx++;
		} while (((cur_loop = BM_vert_step_fan_loop(cur_loop, &cur_edge)) != first_loop) && (cur_loop != NULL));

		//Find the faces for our three points
		{
			int i, j;
			float uvs[3][2];
			BMFace *face_ids[3];
			float rad_dir[3];
			int search_id;

			rad_dir[1] = signf(val_1);
			rad_dir[2] = signf(val_2);

			for ( i = 0; i < edge_count; i++) {
				for ( j = 0; j < 3; j++) {
					float point[3];
					switch(j){
						case 1 :
							copy_v3_v3(point, edge1_mid);
							break;
						case 2 :
							copy_v3_v3(point, edge2_mid);
							break;
						default:
							copy_v3_v3(point, cent);
							break;
					}

					if( point_inside(mat, point, faces[i]) ){
						float point_v2[2];
						float P[3], du[3], dv[3], temp[3];

						mul_v2_m3v3(point_v2, mat, point);

						get_uv_point(faces[i], uvs[j], point_v2, mat);

						face_ids[j] = faces[i];
						if( j == 0 ){
							//Save rad_dir for cent
							openSubdiv_evaluateLimit(m_d->eval, BM_elem_index_get(face_ids[j]), uvs[j][0], uvs[j][1], P, du, dv);

							sub_v3_v3v3(temp, P, C_vert_pos);
							rad_dir[j] = signf(dot_v3v3(rad_plane_no, temp));
						}
					}
				}
			}
			if( rad_dir[0] == rad_dir[1] ){
				search_id = 2;
			} else {
				search_id = 1;
			}

			{
				float search_val, uv_P[2], P[3], du[3], dv[3], temp[3];
				float step = 0.5f;
				float step_len = 0.25f;
				int i, face_index;
				BMFace *orig_face;
				Vert_buf v_buf;
				float prev_val = 0;
				bool changed_val = false;
				/*
				print_v3("cent", cent);
				print_v3("edge1_mid", edge1_mid);
				print_v3("edge2_mid", edge2_mid);
				print_v3("rad_dir", rad_dir);
				print_v2("UV_cent", uvs[0]);
				print_v2("UV_edge1", uvs[1]);
				print_v2("UV_edge2", uvs[2]);
				*/
				if( face_ids[0] == face_ids[search_id] ){
					//We can work in pure uv space
					//printf("UV space\n");
					orig_face = face_ids[0];
					face_index = BM_elem_index_get(face_ids[0]);
					for( i = 0; i < 20; i++ ){
						interp_v2_v2v2( uv_P, uvs[0], uvs[search_id], step);
						openSubdiv_evaluateLimit(m_d->eval, face_index, uv_P[0], uv_P[1], P, du, dv);

						sub_v3_v3v3(temp, P, C_vert_pos);
						search_val = dot_v3v3(rad_plane_no, temp);

						if( fabs(search_val) < 1e-14 ){
							//We got lucky and found the zero crossing!
							printf("got lucky\n");
							break;
						}

						search_val = signf(search_val);

						if( signf(search_val) == rad_dir[0] ){
							step += step_len;
						} else {
							step -= step_len;
						}

						if( !(prev_val == 0) && prev_val != search_val ){
							//Because we are interpolating in coordinate space (and using the limit surface to get new points)
							//we need to be able to go beyond 0 and 1 (so the final step could be 1.2 etc...)
							step_len = step_len/2.0f;
							changed_val = true;
						}
						prev_val = search_val;

						if( i == 9 && !changed_val ){
							//TODO it should not be able to get the wrong edge to search with
							//but it seems like this is the case in some places...
							changed_val = true; //Do not try this again
							i = 0;
							step = 0.5f;
							step_len = 0.25f;
							if(search_id == 2){
								search_id = 1;
							} else {
								search_id = 2;
							}
						}
					}
				} else {
					//Work in coord space
					float cur_p[3], p[3];
					int j;

					//printf("Coord space\n");
					if( search_id == 1 ){
						copy_v3_v3(p, edge1_mid);
					} else {
						copy_v3_v3(p, edge2_mid);
					}

					for( i = 0; i < 20; i++ ){
						interp_v3_v3v3(cur_p, cent, p, step);

						for ( j = 0; j < edge_count; j++) {
							if( point_inside(mat, cur_p, faces[j]) ){
								float point_v2[2];
								mul_v2_m3v3(point_v2, mat, cur_p);

								get_uv_point(faces[j], uv_P, point_v2, mat);

								orig_face = faces[j];
								face_index = BM_elem_index_get(faces[j]);
								openSubdiv_evaluateLimit(m_d->eval, face_index, uv_P[0], uv_P[1], P, du, dv);

								break;
							}
						}

						sub_v3_v3v3(temp, P, C_vert_pos);
						search_val = dot_v3v3(rad_plane_no, temp);

						if( fabs(search_val) < 1e-14 ){
							//We got lucky and found the zero crossing!
							printf("got lucky\n");
							break;
						}

						search_val = signf(search_val);

						if( search_val == rad_dir[0] ){
							step += step_len;
						} else {
							step -= step_len;
						}

						if( !(prev_val == 0) && prev_val != search_val ){
							//Because we are interpolating in coordinate space (and using the limit surface to get new points)
							//we need to be able to go beyond 0 and 1 (so the final step could be 1.2 etc...)
							step_len = step_len/2.0f;
							changed_val = true;
						}
						prev_val = search_val;

						if( i == 9 && !changed_val ){
							//TODO it should not be able to get the wrong edge to search with
							//but it seems like this is the case in some places...
							changed_val = true; //Do not try this again
							i = 0;
							step = 0.5f;
							step_len = 0.25f;
							if(search_id == 2){
								copy_v3_v3(p, edge1_mid);
							} else {
								copy_v3_v3(p, edge2_mid);
							}
						}
					}
				}

				v_buf.orig_edge = NULL;
				v_buf.orig_face = orig_face;
				v_buf.u = uv_P[0];
				v_buf.v = uv_P[1];
				if( poke_and_move(poke_face, P, du, dv, m_d) ){
					BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
				}
			}

		}

	}

}

static void radial_insertion( MeshData *m_d ){

	int vert_i, vert_j;
	BMFace *f;
	BMIter iter_f;
	int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);
	int initial_verts = BM_mesh_elem_count(m_d->bm, BM_VERT);

	for(vert_i = 0; vert_i < m_d->C_verts->count; vert_i++){
		int vert_idx, CC_idx, CC2_idx;
		float co_arr[3][3];
		BMVert *vert_arr[3];
		BMIter iter_v;

		BMVert *vert = BLI_buffer_at(m_d->C_verts, BMVert*, vert_i);

		int face_i;
		int face_count = BM_vert_face_count(vert);
		BMFace *face_arr[face_count];

		/*
		if( BM_elem_index_get(vert) != 1020 ){
			continue;
		}
		*/

		BM_ITER_ELEM_INDEX (f, &iter_f, vert, BM_FACES_OF_VERT, face_i) {
			face_arr[face_i] = f;
		}

		for( face_i = 0; face_i < face_count; face_i++){
			BMVert *cur_vert;
			CC2_idx = -1;
			bool skip_face = false;

			f = face_arr[face_i];

			BM_ITER_ELEM_INDEX (cur_vert, &iter_v, f, BM_VERTS_OF_FACE, vert_idx) {
				if( BM_elem_index_get(cur_vert) > (initial_verts - 1) ){
					// This is a face we already worked on
					skip_face = true;
				}

				for(vert_j = 0; vert_j < m_d->C_verts->count; vert_j++){
					BMVert *vert2 = BLI_buffer_at(m_d->C_verts, BMVert*, vert_j);

					if( cur_vert == vert){
						CC_idx = vert_idx;
						break;
					} else if ( cur_vert == vert2 ){
						CC2_idx = vert_idx;
						break;
					}
				}
				vert_arr[vert_idx] = cur_vert;
				copy_v3_v3(co_arr[vert_idx], cur_vert->co);
			}

			if( skip_face ){
				continue;
			}

			{
				float u_arr[3]; //array for u-coords (v1_u, v2_u ...)
				float v_arr[3];

				float val_1, val_2;
				float temp[3];

				float cam_vec[3], rad_plane_no[3];

				BMFace *orig_face;

				sub_v3_v3v3(cam_vec, m_d->cam_loc, co_arr[CC_idx]);
				cross_v3_v3v3(rad_plane_no, vert_arr[CC_idx]->no, cam_vec);

				{
					//check if the radial plane intersects the vert's opposite edge
					sub_v3_v3v3(temp, co_arr[ mod_i(CC_idx-1, 3) ], co_arr[CC_idx]);

					val_1 = dot_v3v3(rad_plane_no, temp);

					sub_v3_v3v3(temp, co_arr[ mod_i(CC_idx+1, 3) ], co_arr[CC_idx]);

					val_2 = dot_v3v3(rad_plane_no, temp);

					if( signf(val_1) == signf(val_2) ){
						//Do not insert a radial edge here
						continue;
					}
				}

				//get uv coord and orig face
				//TODO rewrite get_orig_face to handle shifted verts.
				//Or just rewrite it in here and leave get_orig_face alone
				orig_face = get_orig_face(orig_verts, vert_arr, u_arr, v_arr, co_arr, m_d);
				if( CC2_idx != -1 ){
					//This face has an CC edge
					//Do the radial planes intersect?
					float rad_plane_no2[3];
					float val2_1, val2_2;

					sub_v3_v3v3(cam_vec, m_d->cam_loc, co_arr[CC2_idx]);
					cross_v3_v3v3(rad_plane_no2, vert_arr[CC2_idx]->no, cam_vec);

					//check if the radial plane intersects the vert's opposite edge
					sub_v3_v3v3(temp, co_arr[ mod_i(CC2_idx-1, 3) ], co_arr[CC2_idx]);

					val2_1 = dot_v3v3(rad_plane_no2, temp);

					sub_v3_v3v3(temp, co_arr[ mod_i(CC2_idx+1, 3) ], co_arr[CC2_idx]);

					val2_2 = dot_v3v3(rad_plane_no2, temp);

					if( signf(val2_1) != signf(val2_2) ){
						//TODO Check if this really works
						//TODO rewrite this to handle shifted verts
						int face_index = BM_elem_index_get(orig_face);
						float uv[2];
						bool b_arr[3];

						b_arr[mod_i(CC_idx-1, 3)] = (val_1 > 0);
						b_arr[mod_i(CC_idx+1, 3)] = (val_2 > 0);
						b_arr[CC_idx] = false;
						printf("Radial intersect!\n");
						if (rad_triangle(m_d->eval, rad_plane_no, rad_plane_no2, co_arr[CC_idx], co_arr[CC2_idx],
									co_arr, b_arr, u_arr, v_arr, face_index, uv))
						{
							float P[3], du[3], dv[3];
							Vert_buf v_buf;
							v_buf.orig_edge = NULL;
							v_buf.orig_face = orig_face;
							v_buf.u = uv[0];
							v_buf.v = uv[1];

							openSubdiv_evaluateLimit(m_d->eval, face_index, uv[0], uv[1], P, du, dv);
							if( poke_and_move(f, P, du, dv, m_d) ){
								BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
							}
							printf("Found radi point!\n");
							continue;
						}
					}
				}

				{
					int i;
					bool diff_faces = false;
					BMFace *faces[3];
					//Check if the triangle has been shifted so we can't use the original face for UV coords
					for(i = 0; i < 3; i++){
						Vert_buf* vert = get_shift_vert( vert_arr[i], m_d );
						if( vert != NULL ){
							//This vert has been shifted
							if(vert->orig_face != orig_face){
								//We have shifted this vert into an other face.
								//Now we have to interpolate on multiple faces
								diff_faces = true;
							}
							faces[i] = vert->orig_face;
							//Update to new uv pos
							u_arr[i] = vert->u;
							v_arr[i] = vert->v;
						} else {
							//Check if edge verts doesn't belong to orig_face
							int v_idx = BM_elem_index_get(vert_arr[i]);
							if( (v_idx + 1) > orig_verts){
								Vert_buf v_buf = BLI_buffer_at(m_d->new_vert_buffer, Vert_buf, v_idx - orig_verts);
								if( v_buf.orig_edge != NULL ){
									BMIter iter_f;
									BMFace *face;
									bool found_face = false;

									//TODO try to find a edge face that is next to orig_face if none of them is orig_face
									BM_ITER_ELEM (face, &iter_f, v_buf.orig_edge, BM_FACES_OF_EDGE){
										if(face == orig_face){
											found_face = true;
											break;
										}
									}

									if(!found_face){
										faces[i] = v_buf.orig_face;
										u_arr[i] = v_buf.u;
										v_arr[i] = v_buf.v;
										diff_faces = true;
										continue;
									}
								}
							}
							faces[i] = orig_face;
						}
					}

					if( true ){
						//This search will spawn multiple faces, we must use coordinate space to do this.
						float cent[3];
						float edge1_mid[3];
						float edge2_mid[3];
						BMFace *face_idx1,*face_idx2;
						interp_v3_v3v3(edge1_mid, co_arr[CC_idx], co_arr[ mod_i(CC_idx-1, 3) ], 0.5f);
						interp_v3_v3v3(edge2_mid, co_arr[CC_idx], co_arr[ mod_i(CC_idx+1, 3) ], 0.5f);
						cent_tri_v3(cent, co_arr[0], co_arr[1], co_arr[2]);

						printf("Diff faces\n");
						mult_radi_search(faces, cent, edge1_mid, edge2_mid, val_1, val_2, rad_plane_no, co_arr[CC_idx], f, m_d);
						continue;
					}

				}

				{
					//get center uv coords
					float cent_uv[2] = {(u_arr[0] + u_arr[1] + u_arr[2]) / 3.0f, (v_arr[0] + v_arr[1] + v_arr[2]) / 3.0f};
					//get edge uv midpoints
					float edge1_mid[2] = { 0.5f * (u_arr[CC_idx] + u_arr[ mod_i(CC_idx-1, 3) ]),
										0.5f * (v_arr[CC_idx] + v_arr[ mod_i(CC_idx-1, 3) ]) };
					float edge2_mid[2] = { 0.5f * (u_arr[CC_idx] + u_arr[ mod_i(CC_idx+1, 3) ]),
										0.5f * (v_arr[CC_idx] + v_arr[ mod_i(CC_idx+1, 3) ]) };

					int	face_index = BM_elem_index_get(orig_face);
					float du[3], dv[3], P[3];
					float search_val;
					//Begin search
					openSubdiv_evaluateLimit(m_d->eval, face_index, cent_uv[0], cent_uv[1], P, du, dv);

					sub_v3_v3v3(temp, P, co_arr[CC_idx]);
					search_val = dot_v3v3(rad_plane_no, temp);
					//Did we get lucky?
					if( fabs(search_val) < 1e-14 ){
						//Insert new vertex
						Vert_buf v_buf;
						v_buf.orig_edge = NULL;
						v_buf.orig_face = orig_face;
						v_buf.u = cent_uv[0];
						v_buf.v = cent_uv[1];
						if( poke_and_move(f, P, du, dv, m_d) ){
							BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
						}
						printf("got lucky\n");
						continue;
					}

					{
						int i;
						float uv_P[2];
						float step = 0.5f;
						float step_len = 0.25f;
						if( signf(search_val) != signf(val_1) ){
							for( i = 0; i < 10; i++){
								interp_v2_v2v2( uv_P, edge1_mid, cent_uv, step);
								openSubdiv_evaluateLimit(m_d->eval, face_index, uv_P[0], uv_P[1], P, du, dv);

								sub_v3_v3v3(temp, P, co_arr[CC_idx]);
								search_val = dot_v3v3(rad_plane_no, temp);

								if( fabs(search_val) < 1e-14 ){
									//We got lucky and found the zero crossing!
									printf("got lucky\n");
									break;
								}

								if( signf(search_val) == signf(val_1) ){
									step += step_len;
								} else {
									step -= step_len;
								}
								step_len = step_len/2.0f;
							}
							{
								Vert_buf v_buf;
								v_buf.orig_edge = NULL;
								v_buf.orig_face = orig_face;
								v_buf.u = cent_uv[0];
								v_buf.v = cent_uv[1];
								if( poke_and_move(f, P, du, dv, m_d) ){
									BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
								}
							}
						} else {
							for( i = 0; i < 10; i++){
								interp_v2_v2v2( uv_P, edge2_mid, cent_uv, step);
								openSubdiv_evaluateLimit(m_d->eval, face_index, uv_P[0], uv_P[1], P, du, dv);

								sub_v3_v3v3(temp, P, co_arr[CC_idx]);
								search_val = dot_v3v3(rad_plane_no, temp);

								if( fabs(search_val) < 1e-14 ){
									//We got lucky and found the zero crossing!
									printf("got lucky\n");
									break;
								}

								if( signf(search_val) == signf(val_2) ){
									step += step_len;
								} else {
									step -= step_len;
								}
								step_len = step_len/2.0f;
							}
							{
								Vert_buf v_buf;
								v_buf.orig_edge = NULL;
								v_buf.orig_face = orig_face;
								v_buf.u = cent_uv[0];
								v_buf.v = cent_uv[1];
								if( poke_and_move(f, P, du, dv, m_d) ){
									BLI_buffer_append(m_d->new_vert_buffer, Vert_buf, v_buf);
								}
							}
						}
					}
				}
			}
		}
	}
}

static bool radial_C_vert(BMVert *v, int * const radi_start_idx, BLI_Buffer *C_verts){


	if( !(BM_elem_index_get(v) < *radi_start_idx) ){
		//This is a radial vert
		return true;
	}

	if( is_C_vert( v, C_verts) ){
		return true;
	}

	return false;
}

static void radial_flip( MeshData *m_d ){

	bool done = false;

	int iter = 0;

	while( !done ){

	int flips = 0;
	int vert_i;
	for(vert_i = 0; vert_i < m_d->C_verts->count; vert_i++){
		BMEdge *e;
		BMIter iter_e;
		int edge_idx;
		int edge_i;
		BMVert *vert = BLI_buffer_at(m_d->C_verts, BMVert*, vert_i);
		int edge_count = BM_vert_edge_count(vert);
		BMEdge *edge_arr[edge_count];

		BM_ITER_ELEM_INDEX (e, &iter_e, vert, BM_EDGES_OF_VERT, edge_idx) {
			edge_arr[edge_idx] = e;
		}

		for( edge_i = 0; edge_i < edge_count; edge_i++){
			BMVert *edge_vert;
			e = edge_arr[edge_i];

			if(e->v1 != vert){
				edge_vert = e->v1;
			} else {
				edge_vert = e->v2;
			}

			if( radial_C_vert( edge_vert, &(m_d->radi_start_idx), m_d->C_verts ) ){
				//This is a radial or CC edge, do not try to flip it.
				continue;
			}

			//See if it's possible to rotate the edge at all
			// IE check for mesh border edges etc
			if( !BM_edge_rotate_check(e) ){
				//Not possible, bail
				printf("Couldn't rotate edge!\n");
				continue;
			}

			{
				//Check if we can just do a simple rotation that doesn't create any bad faces
				BMLoop *l1, *l2;
				BM_edge_calc_rotate(e, true, &l1, &l2);

				// new_v1 = l1->v;
				// new_v2 = l2->v;

				if( BM_elem_index_get(l1->v) < m_d->radi_start_idx &&
					BM_elem_index_get(l2->v) < m_d->radi_start_idx ){
					//The flip will not increase the number of standard radial triangles
					//Do not flip it
					continue;
				}
				//Check if the flip creates any folds
				{
					float mat[3][3];
					float mat_coords[3][2], mat_new_pos[2];

					axis_dominant_v3_to_m3(mat, edge_vert->no);
					mul_v2_m3v3(mat_new_pos, mat, edge_vert->co);

					mul_v2_m3v3(mat_coords[0], mat, vert->co);
					mul_v2_m3v3(mat_coords[1], mat, l1->v->co);
					mul_v2_m3v3(mat_coords[2], mat, l2->v->co);

					if( !isect_point_tri_v2(mat_new_pos, mat_coords[0], mat_coords[1], mat_coords[2]) ){

						//We can simply rotate it!
						BM_edge_rotate(m_d->bm, e, true, 0);
						flips++;
						continue;
					}
				}

				printf("Try to dissolve vert!\n");

				{
					int edge_count = BM_vert_edge_count(edge_vert);
					BMEdge *cur_e = e;
					BMEdge *edge_arr[edge_count];
					int edge_idx = 0;

					BMEdge *rad1_edge = BM_edge_exists(l1->v, edge_vert);
					BMEdge *rad2_edge = BM_edge_exists(l2->v, edge_vert);

					BMLoop *first_loop = BM_face_vert_share_loop( cur_e->l->f, edge_vert);
					BMLoop *cur_loop = first_loop;

					if( edge_count - 3 < 1 ){
						continue;
					}

					if( rad1_edge == NULL || rad2_edge == NULL){
						printf("Couldn't find the edges connected to the radial vert! (dissolve vert)\n");
						continue;
					}

					if( cur_loop == NULL ){
						printf("Couldn't find face loop!\n");
					}

					while (((cur_loop = BM_vert_step_fan_loop(cur_loop, &cur_e)) != first_loop) && (cur_loop != NULL)) {
						if(cur_e == e || cur_e == rad1_edge || cur_e == rad2_edge){
							continue;
						}

						edge_arr[edge_idx] = cur_e;
						edge_idx++;

					}

					if(cur_loop == NULL){
						continue;
					}

					for( edge_idx = 0; edge_idx < edge_count - 3; edge_idx++){
						BMLoop *loop1, *loop2;
						BM_edge_calc_rotate(edge_arr[edge_idx], true, &loop1, &loop2);

						if( BM_edge_rotate_check_degenerate(edge_arr[edge_idx], loop1, loop2) ){
							BM_edge_rotate(m_d->bm, edge_arr[edge_idx], true, 0);
						} else {
							//Try to rotate from the other side instead
							printf("Try from other side!\n");
							break;
						}
					}

					if( edge_idx != edge_count - 3 ){
						int op_idx = edge_count -4;
						bool failed_rotate = false;

						for(; edge_idx <= op_idx; op_idx--){

							BMLoop *loop1, *loop2;
							BM_edge_calc_rotate(edge_arr[op_idx], true, &loop1, &loop2);

							if( edge_idx == op_idx ){
								//This is the last edge that is going to be rotated
								// check_degenerate is too strict in this case (in most cases the face area will be near zero)
								//TODO check for folds
								if(	BM_edge_rotate(m_d->bm, edge_arr[op_idx], true, 0) != NULL ){
									continue;
								} else {
									failed_rotate = true;
									break;
								}
							}

							if( BM_edge_rotate_check_degenerate(edge_arr[op_idx], loop1, loop2) ){
								BM_edge_rotate(m_d->bm, edge_arr[op_idx], true, 0);
							} else {
								failed_rotate = true;
								break;
							}
						}

						if( failed_rotate ){
							printf("Failed to flip and delete in radial edge flip!\n");
							continue;
						}

					}
				}

				if( BM_disk_dissolve(m_d->bm, edge_vert) ){
					//Count down m_d->radi_start_idx because we have fewer verts after the dissolve;
					flips++;
					m_d->radi_start_idx -= 1;
					printf("Dissolved vert!\n");
				}

			}

		}
	}

	if(flips == 0 || iter > 5){
		done = true;
	}
	iter++;
	}
}

static void optimization( MeshData *m_d ){

	BLI_buffer_declare_static(IncoFace, inco_faces, BLI_BUFFER_NOP, 32);
	//Find and save all inconsistent faces before we begin trying to optimize the mesh
	{
		int i;
		BMVert *vert;
		BMIter iter_v;
		//TODO is there anyway to begin at m_d->radi_start_idx?
		BM_ITER_MESH_INDEX (vert, &iter_v, m_d->bm, BM_VERTS_OF_MESH, i) {

			if( i < m_d->radi_start_idx ){
				continue;
			}

			{
				BMFace *face;
				BMIter iter_f;
				bool b_f = calc_if_B_nor(m_d->cam_loc, vert->co, vert->no);
				float P[3];

				BM_ITER_ELEM (face, &iter_f, vert, BM_FACES_OF_VERT) {
					//TODO mark inconsistent faces in an other way
					// and only check each face once
					if(face->mat_nr == 4){
						//Already added this face to inco_faces
						continue;
					}

					{
						// This shouldn't be needed, but in case we manage to have inconsistent
						// faces that borders our contour line, don't mark it for adjustment.
						BMVert *v;
						BMIter iter;
						bool found_c_vert = false;

						BM_ITER_ELEM (v, &iter, face, BM_VERTS_OF_FACE) {
							if( is_C_vert(v, m_d->C_verts) ) {
								found_c_vert = true;
								break;
							}
						}

						if( found_c_vert ) {
							continue;
						}

					}
					BM_face_calc_center_mean(face, P);

					if( b_f != calc_if_B_nor(m_d->cam_loc, P, face->no) ){
						IncoFace inface;
						inface.face = face;
						inface.back_f = b_f;
						face->mat_nr = 4;
						BLI_buffer_append(&inco_faces, IncoFace, inface);
					}

				}
			}

		}
	}
	// 1. Radial edge extension
	{
		// TODO
	}

	// 2. Edge flipping
	{
		int face_i;

		for(face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);

			BMEdge *edge;
			BMIter iter_e;

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}

			BM_ITER_ELEM (edge, &iter_e, inface->face, BM_EDGES_OF_FACE) {

				if( !BM_edge_rotate_check(edge) ){
					continue;
				}

				if( BM_elem_index_get(edge->v1) < m_d->radi_start_idx || BM_elem_index_get(edge->v2) < m_d->radi_start_idx ){
					//This is not a radial triangle edge, see if we can flip it
					BMLoop *l1, *l2;
					BM_edge_calc_rotate(edge, true, &l1, &l2);

					if( !BM_edge_rotate_check_degenerate(edge, l1, l2) ){
						continue;
					}

					if( !BM_edge_rotate_check_beauty(edge, l1, l2) ){
						continue;
					}

					{
						float vec1[3], vec2[3], P[3], no[3];
						BMVert *v1, *v2;

						BM_edge_ordered_verts(edge, &v1, &v2);

						//TODO perhaps use normal_tri_v3 instead for normal calc

						sub_v3_v3v3(vec1, v1->co, l1->v->co);
						sub_v3_v3v3(vec2, v1->co, l2->v->co);

						cross_v3_v3v3(no, vec1, vec2);
						normalize_v3(no);

						//Calc center mean of new face
						zero_v3(P);
						add_v3_v3( P, v1->co );
						add_v3_v3( P, l1->v->co );
						add_v3_v3( P, l2->v->co );

						mul_v3_fl( P, 1.0f / 3.0f );

						if( inface->back_f != calc_if_B_nor(m_d->cam_loc, P, no) ){
							//This is not a good flip!
							printf("Opti flip, first face not good\n");
							continue;
						}

						sub_v3_v3v3(vec1, v2->co, l1->v->co);
						sub_v3_v3v3(vec2, v2->co, l2->v->co);

						cross_v3_v3v3(no, vec2, vec1);
						normalize_v3(no);

						//Calc center mean of new face
						zero_v3(P);
						add_v3_v3( P, v2->co );
						add_v3_v3( P, l1->v->co );
						add_v3_v3( P, l2->v->co );

						mul_v3_fl( P, 1.0f / 3.0f );

						if( inface->back_f != calc_if_B_nor(m_d->cam_loc, P, no) ){
							//This is not a good flip!
							printf("Opti flip, second face not good\n");
							continue;
						}

						printf("Opti filped an edge!\n");

						//TODO remove this or only when debug
						inface->face->mat_nr = 0;

						BM_edge_rotate(m_d->bm, edge, true, 0);
						inface->face = NULL;
						//Done with this face
						break;
					}
				}

			}

		}
	}
	// 2.a (Not in the paper) Vertex dissolve
	{
		int face_i;

		for(face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);

			BMVert *vert;
			BMIter iter_v;

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}

			BM_ITER_ELEM (vert, &iter_v, inface->face, BM_VERTS_OF_FACE) {
				if( BM_elem_index_get(vert) < m_d->radi_start_idx ){
					//Not a radial vert, see if we can dissolve it to improve the consistency

					if( BM_vert_edge_count(vert) == 3 && BM_vert_face_count(vert) == 3 ){
						float vec1[3], vec2[3], P[3], no[3];
						BMLoop *l1, *l2;
						BMVert *v1, *v2;
						BM_edge_calc_rotate(vert->e, true, &l1, &l2);

						BM_edge_ordered_verts(vert->e, &v1, &v2);

						if( vert == v1 ){
							sub_v3_v3v3(vec1, v2->co, l1->v->co);
							sub_v3_v3v3(vec2, v2->co, l2->v->co);

							cross_v3_v3v3(no, vec2, vec1);
							normalize_v3(no);

							//Calc center mean of new face
							zero_v3(P);
							add_v3_v3( P, v2->co );
							add_v3_v3( P, l1->v->co );
							add_v3_v3( P, l2->v->co );

							mul_v3_fl( P, 1.0f / 3.0f );
						} else {
							sub_v3_v3v3(vec1, v1->co, l1->v->co);
							sub_v3_v3v3(vec2, v1->co, l2->v->co);

							cross_v3_v3v3(no, vec1, vec2);
							normalize_v3(no);

							//Calc center mean of new face
							zero_v3(P);
							add_v3_v3( P, v1->co );
							add_v3_v3( P, l1->v->co );
							add_v3_v3( P, l2->v->co );

							mul_v3_fl( P, 1.0f / 3.0f );
						}

						if(inface->back_f == calc_if_B_nor(m_d->cam_loc, P, no)){
							printf("Opti dissolve\n");
							//TODO remove this or only when debug
							inface->face->mat_nr = 0;

							if(!BM_disk_dissolve(m_d->bm, vert)){
								printf("Failed to opti dissolve\n");
								continue;
							}
							//Count down radi_start_idx because we have fewer verts after the dissolve;
							m_d->radi_start_idx -= 1;

							inface->face = NULL;
							//Done with this face
							break;
						}

					}

				}
			}
		}
	}

	// 2.b (Not in the paper) Smooth vertex position
	// TODO perhaps move this to before wiggling in normal direction (IE after step 4)
	{
		int face_i;

		for(face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);

			BMVert *vert;
			BMIter iter_v;

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}

			BM_ITER_ELEM (vert, &iter_v, inface->face, BM_VERTS_OF_FACE) {
				if( BM_elem_index_get(vert) < m_d->radi_start_idx ){
					//not a radial vert, try to smooth the vertex pos and see if the consistency improves

					float old_pos[3], co[3], co2[3];
					int i = 0;
					bool done = true;

					BMEdge *edge;
					BMIter iter_e;

					zero_v3(co);
					copy_v3_v3(old_pos, vert->co);


					BM_ITER_ELEM (edge, &iter_e, vert, BM_EDGES_OF_VERT) {
						copy_v3_v3(co2, BM_edge_other_vert(edge, vert)->co);
						add_v3_v3v3(co, co, co2);
						i += 1;
					}

					mul_v3_fl(co, 1.0f / (float)i);
					mid_v3_v3v3(co, co, vert->co);

					copy_v3_v3(vert->co, co);

					{
						BMFace *face;
						BMIter iter_f;

						BM_ITER_ELEM (face, &iter_f, vert, BM_FACES_OF_VERT) {
							float no[3];
							float P[3];
							BM_face_calc_normal(face, no);
							BM_face_calc_center_mean(face, P);

							if( inface->back_f != calc_if_B_nor(m_d->cam_loc, P, no) ){
								//Bad vertex move
								printf("Bad vert smooth\n");
								//Move the vertex back to it's original position

								copy_v3_v3(vert->co, old_pos);
								done = false;
								break;
							}

						}

						if( done ){
							//Good vert smooth
							inface->face->mat_nr = 0;
							inface->face = NULL;
							break;
						}
					}

				}
			}
		}
	}

	// 3. Vertex wiggling in paramter space
	/*{
		int face_i;

		for(face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);
			int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);

			BMVert *vert;
			BMIter iter_v;

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}

			BM_ITER_ELEM (vert, &iter_v, inface->face, BM_VERTS_OF_FACE) {
				int idx = BM_elem_index_get(vert);
				if( idx < orig_verts ){
					// This vert exists in the original mesh
					BMVert *orig_v = BM_vert_at_index_find(m_d->bm_orig, idx);
					int edge_count = BM_vert_edge_count(vert);

					// We are going to search for a new position using a spiral
					float end_radius = 0;
					float cur_radian = 0;
					float max_turns = 4; // How many turnings the spiral should have

					bool done = false;

					BMEdge *e;
					BMIter iter;
					BM_ITER_ELEM (e, &iter, vert, BM_EDGES_OF_VERT) {
						end_radius += BM_edge_calc_length(e);
					}

					end_radius = end_radius / (float)edge_count;
					end_radius = end_radius / 2.0f;

					{
						// radius = b * theta
						float b = (end_radius / max_turns) * (1.0f / (2.0f * M_PI));
						float step_size = (max_turns * 2.0f * M_PI) / 100.0f;
						int i;

						BMFace *f;
						BMIter iter_f;

						float mat[3][3];
						float center_v2[2];
						float old_pos[3];

						copy_v3_v3( old_pos, vert->co );

						axis_dominant_v3_to_m3(mat, vert->no);
						mul_v2_m3v3(center_v2, mat, vert->co);

						for( i=0; i<100; i++ ){
							float theta = step_size * (float)i;
							float r = b * theta;
							float pos_v2[2] = { center_v2[0], center_v2[1] };

							bool found_point = true;

							pos_v2[0] += r * cosf(theta);
							pos_v2[1] += r * sinf(theta);

							// TODO check if the new point lies inside any of the new mesh faces

							BM_ITER_ELEM (f, &iter_f, orig_v, BM_FACES_OF_VERT) {
								if( point_inside_v2( mat, pos_v2, f ) ){
									float P[3], du[3], dv[3];
									float uv_P[2];

									get_uv_point( f, uv_P, pos_v2, mat );
									openSubdiv_evaluateLimit(m_d->eval, BM_elem_index_get(f), uv_P[0], uv_P[1], P, du, dv);

									copy_v3_v3(vert->co, P);
									//No need to iterate over the remaining faces
									break;
								}
							}

							BM_ITER_ELEM (f, &iter_f, vert, BM_FACES_OF_VERT) {
								float no[3];
								float P[3];
								BM_face_calc_normal(f, no);
								BM_face_calc_center_mean(f, P);

								if( inface->back_f != calc_if_B_nor(m_d->cam_loc, P, no) ){
									found_point = false;
									break;
								}
							}

							if( found_point ){
								done = true;
								break;
							}
						}

						if( done ){
							inface->face->mat_nr = 0;
							inface->face = NULL;
							printf("Vertex wiggle\n");
							break;
						} else {
							printf("Bad Vertex wiggle\n");
							copy_v3_v3(vert->co, old_pos);
						}
					}


				}
			}
		}
	}*/

	// 4. Edge Splitting
	/*{
		// TODO
		int face_i;

		for(face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);
			int orig_verts = BM_mesh_elem_count(m_d->bm_orig, BM_VERT);

			BMEdge *edge;
			BMIter iter_v;

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}

			BM_ITER_ELEM (edge, &iter_v, inface->face, BM_EDGES_OF_FACE) {
				BMVert *orig_v = NULL;
				BMFace *f;
				BMIter iter_f;

				// Storage for vertex pos
				// This is used later to see if the edge split will be ok or not
				float store[4][2][3];
				int j = 0;

				if( BM_elem_index_get(edge->v1) < orig_verts ){
					int idx = BM_elem_index_get(edge->v1);
					orig_v = BM_vert_at_index_find(m_d->bm_orig, idx);
				} else if ( BM_elem_index_get(edge->v2) < orig_verts ){
					int idx = BM_elem_index_get(edge->v2);
					orig_v = BM_vert_at_index_find(m_d->bm_orig, idx);
				}

				if( orig_v == NULL ){
					continue;
				}

				BM_ITER_ELEM (f, &iter_f, edge, BM_FACES_OF_EDGE) {
					int i;
					BMLoop *l;
					l = BM_FACE_FIRST_LOOP(f);
					for( i = 0; i < 3; i++ ){
						if( l->e != edge ) {
							copy_v3_v3(store[j][0], l->v->co);
							copy_v3_v3(store[j][1], (l->next)->v->co);
							j++;
						}
						l = l->next;
					}
				}

				{
					//Do 10 samples but don't check end and start point
					float step = 1.0f/11.0f;
					float step_arr[] = { step*5.0f, step*6.0f, step*4.0f, step*7.0f, step*3.0f,
						step*8.0f, step*2.0f, step*9.0f, step*1.0f, step*10.0f };


					float mat[3][3];
					float start[2], end[2];
					int i;
					bool done = false;

					axis_dominant_v3_to_m3(mat, orig_v->no);
					mul_v2_m3v3(start, mat, edge->v1->co);
					mul_v2_m3v3(end, mat, edge->v2->co);
					for( i = 0; i < 10; i++ ){
						float cur_v2[2];
						float P[3], du[3], dv[3];
						bool found_point = false;

						interp_v2_v2v2(cur_v2, start, end, step_arr[i]);

						BM_ITER_ELEM (f, &iter_f, orig_v, BM_FACES_OF_VERT) {
							if( point_inside_v2( mat, cur_v2, f ) ){
								float uv_P[2];

								get_uv_point( f, uv_P, cur_v2, mat );
								openSubdiv_evaluateLimit(m_d->eval, BM_elem_index_get(f), uv_P[0], uv_P[1], P, du, dv);

								found_point = true;
								//No need to iterate over the remaining faces
								break;
							}
						}

						if( found_point ) {
							float p_cent[3];
							float no[3];
							done = true;

							for( j = 0; j < 4; j++ ){
								zero_v3(p_cent);
								add_v3_v3(p_cent, store[j][0]);
								add_v3_v3(p_cent, store[j][1]);
								add_v3_v3(p_cent, P);

								mul_v3_fl(p_cent, 1.0f / 3.0f);

								normal_tri_v3(no, store[j][0], store[j][1], P);
								if( inface->back_f != calc_if_B_nor(m_d->cam_loc, p_cent, no) ){
									done = false;
									break;
								}
							}

							if( done ) {
								inface->face->mat_nr = 0;
								inface->face = NULL;
								print_v3("P", P);
								//TODO this will fail because we have dissolved stuff before...
								split_edge_and_move_vert(m_d->bm, edge, P, du, dv);

								break;
							}
						}

					}
					if( done ) {
						break;
					}

				}
			}
		}
	}*/

	// 5. Vertex wiggling in normal direction
	{
		int face_i;

		for(face_i = 0; face_i < inco_faces.count; face_i++){
			IncoFace *inface = &BLI_buffer_at(&inco_faces, IncoFace, face_i);

			BMVert *vert;
			BMIter iter_v;

			if( inface->face == NULL ){
				//Already fixed this edge
				continue;
			}

			BM_ITER_ELEM (vert, &iter_v, inface->face, BM_VERTS_OF_FACE) {
				if( BM_elem_index_get(vert) < m_d->radi_start_idx ){
					BMEdge *edge;
					BMIter iter_e;
					float old_pos[3];
					float len = 0.0f;
					int i = 0;
					bool done = false;

					copy_v3_v3(old_pos, vert->co);

					BM_ITER_ELEM (edge, &iter_e, vert, BM_EDGES_OF_VERT) {
						len += BM_edge_calc_length(edge);
						i++;
					}

					len = len / (float)(i*2);

					for( i = 0; i < 100; i++ ){
						BMFace *face;
						BMIter iter_f;
						bool found_point = true;
						float co[3];
						float cur_len;

						copy_v3_v3(co, vert->no);

						if( i % 2 ){
							cur_len = len * (float)((i+1)/2) / 50.0f;
						} else {
							cur_len = len * (float)((i+2)/2) / -50.0f;
						}

						mul_v3_fl(co, cur_len);

						add_v3_v3v3(vert->co, old_pos, co);

						BM_ITER_ELEM (face, &iter_f, vert, BM_FACES_OF_VERT) {
							float no[3];
							float P[3];
							BM_face_calc_normal(face, no);
							BM_face_calc_center_mean(face, P);

							if( inface->back_f != calc_if_B_nor(m_d->cam_loc, P, no) ){
								found_point = false;
								break;
							}
						}

						if( found_point ){
							done = true;
							break;
						}

					}

					if( done ){
						inface->face->mat_nr = 0;
						inface->face = NULL;
						printf("Opti normal wiggle\n");
						break;
					} else {
						copy_v3_v3(vert->co, old_pos);
					}
				}
			}
		}
	}
	//Cleanup
	BLI_buffer_free(&inco_faces);
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
		if( f->mat_nr == 4 ){
			continue;
		}
		BM_face_calc_center_mean(f, P);
		if( calc_if_B_nor(cam_loc, P, f->no) ){
			f->mat_nr = 1;
		} else {
			f->mat_nr = 0;
		}
	}
}

static void debug_colorize_radi( MeshData *m_d ){
	BMIter iter, iter_e;
	BMFace *f;
	BMEdge *e;
	float P[3];

	int vert_i;
	for(vert_i = 0; vert_i < m_d->C_verts->count; vert_i++){
		BMVert *vert = BLI_buffer_at(m_d->C_verts, BMVert*, vert_i);

		BM_ITER_ELEM (e, &iter_e, vert, BM_EDGES_OF_VERT) {
			BMVert *edge_vert;

			if(e->v1 != vert){
				edge_vert = e->v1;
			} else {
				edge_vert = e->v2;
			}

			if( !(BM_elem_index_get(edge_vert) < m_d->radi_start_idx) ){
				//This is a radial/CC edge vert.
				BM_ITER_ELEM (f, &iter, e, BM_FACES_OF_EDGE){
					if( f->mat_nr == 4 ){
						continue;
					}
					BM_face_calc_center_mean(f, P);
					if( calc_if_B_nor(m_d->cam_loc, P, f->no) ){
						f->mat_nr = 3;
					} else {
						f->mat_nr = 2;
					}
				}
			}
		}
	}
}

/* bmesh only function */
static DerivedMesh *mybmesh_do(DerivedMesh *dm, MyBMeshModifierData *mmd, float cam_loc[3])
{

	DerivedMesh *result;
	BMesh *bm_orig, *bm;
	bool quad_mesh = true;

	//CCGSubSurf *ss;
	struct OpenSubdiv_EvaluatorDescr *osd_eval;

	bm = DM_to_bmesh(dm, true);

	//TODO use this to check if we need to subdivide the mesh to get a quad mesh.
	{
		BMIter iter;
		BMFace *f, *f_next;
		int sides = 4;
		/* use the mutable iterator so we can remove data as its looped over */
		BM_ITER_MESH_MUTABLE (f, f_next, &iter, bm, BM_FACES_OF_MESH) {
			if (f->len != sides) {
				quad_mesh = false;
				break;
			}
		}
	}

	if (!quad_mesh){
		if( mmd->camera_ob != NULL ){
			debug_colorize(bm, cam_loc);
		}

		result = CDDM_from_bmesh(bm, true);

		BM_mesh_free(bm);

		result->dirty |= DM_DIRTY_NORMALS;

		return result;
	}

	//Keep a copy of the original mesh
	bm_orig = DM_to_bmesh(dm, true);

	osd_eval = create_osd_eval(bm);

	// (6.1) Initialization
	verts_to_limit(bm, osd_eval);

	if (mmd->flag & MOD_MYBMESH_TRIANG) {
		//TODO check if shortest diagonal is better
		BM_mesh_triangulate(bm, MOD_TRIANGULATE_QUAD_FIXED, MOD_TRIANGULATE_NGON_BEAUTY, false, NULL, NULL);
	}

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
		BLI_buffer_declare_static(Vert_buf, shifted_verts, BLI_BUFFER_NOP, 32);
		BLI_buffer_declare_static(Cusp, cusp_edges, BLI_BUFFER_NOP, 32);
		BLI_buffer_declare_static(BMVert*, C_verts, BLI_BUFFER_NOP, 32);

		MeshData mesh_data;

		mesh_data.bm = bm;
		mesh_data.bm_orig = bm_orig;

		copy_v3_v3( mesh_data.cam_loc, cam_loc );

		mesh_data.new_vert_buffer = &new_vert_buffer;
		mesh_data.shifted_verts = &shifted_verts;
		mesh_data.cusp_edges = &cusp_edges;
		mesh_data.C_verts = &C_verts;
		mesh_data.is_cusp = false;
		mesh_data.eval = osd_eval;

		if (mmd->flag & MOD_MYBMESH_FF_SPLIT) {
			split_BB_FF_edges(&mesh_data);
		}
		// (6.2) Contour Insertion

		//TODO implement vertex shift (as an alternative to edge split)

		if (mmd->flag & MOD_MYBMESH_CUSP_D) {
			cusp_detection(&mesh_data);
		}

		if (mmd->flag & MOD_MYBMESH_FB_SPLIT) {
			contour_insertion(&mesh_data);
		}
		if (mmd->flag & MOD_MYBMESH_CUSP_I) {
			cusp_insertion(&mesh_data);
		}

		// (6.3) Radialization

		mesh_data.radi_start_idx = BM_mesh_elem_count(bm, BM_VERT);

		if (mmd->flag & MOD_MYBMESH_RAD_I){
			radial_insertion(&mesh_data);
		}

		if (mmd->flag & MOD_MYBMESH_RAD_FLIP){
			radial_flip(&mesh_data);
		}

		//Recalculate normals
		BM_mesh_normals_update(bm);

		// (6.4) Optimization
		if (mmd->flag & MOD_MYBMESH_OPTI){
			optimization(&mesh_data);
			//Recalculate normals
			BM_mesh_normals_update(bm);
		}

		debug_colorize(bm, cam_loc);
		debug_colorize_radi(&mesh_data);
		BLI_buffer_free(&new_vert_buffer);
		BLI_buffer_free(&shifted_verts);
		BLI_buffer_free(&cusp_edges);
		BLI_buffer_free(&C_verts);
	}
	result = CDDM_from_bmesh(bm, true);

	BM_mesh_free(bm);
	BM_mesh_free(bm_orig);

	openSubdiv_deleteEvaluatorDescr(osd_eval);

	result->dirty |= DM_DIRTY_NORMALS;

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
		/*
		printf("Cam loc:\n");
		printf("1: %f\n", cam_loc[0]);
		printf("2: %f\n", cam_loc[1]);
		printf("3: %f\n", cam_loc[2]);
		*/
		//convert camera origin from world coord to the modifier obj local coords
		mul_m4_v3(ob->obmat, cam_loc);
		/*
		printf("Cam loc 2:\n");
		printf("1: %f\n", cam_loc[0]);
		printf("2: %f\n", cam_loc[1]);
		printf("3: %f\n", cam_loc[2]);
		*/
	}

	if (!(result = mybmesh_do(dm, mmd, cam_loc))) {
		return dm;
	}

	return result;
}

static void foreachObjectLink(ModifierData *md, Object *ob,
							void (*walk)(void *userData, Object *ob, Object **obpoin),
							void *userData)
{
	MyBMeshModifierData *mmd = (MyBMeshModifierData *)md;

	walk(userData, ob, &mmd->camera_ob);
}

static void updateDepgraph(ModifierData *md, DagForest *forest,
						struct Scene *UNUSED(scene),
						Object *UNUSED(ob),
						DagNode *obNode)
{
	MyBMeshModifierData *mmd = (MyBMeshModifierData *)md;

	if (mmd->camera_ob) {
		DagNode *latNode = dag_get_node(forest, mmd->camera_ob);

		dag_add_relation(forest, latNode, obNode,
				DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "MyBmesh Modifier");
	}
}

static bool dependsOnNormals(ModifierData *UNUSED(md))
{
	//TODO it does depend on normals. return true here?
	return false;
}

ModifierTypeInfo modifierType_MyBMesh = {
	/* name */              "MyBMesh",
	/* structName */        "MyBMeshModifierData",
	/* structSize */        sizeof(MyBMeshModifierData),
	/* type */              eModifierTypeType_Constructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	/* flags cont */        eModifierTypeFlag_SupportsEditmode |
	/* flags cont */        eModifierTypeFlag_EnableInEditmode,

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
	/* updateDepgraph */    updateDepgraph,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */  dependsOnNormals,
	/* foreachObjectLink */ foreachObjectLink,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};
