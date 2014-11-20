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

		//TODO can I just check each very once?
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

	split_edge_and_move_vert(bm, e, P, du, dv, uv_P[0], uv_P[1]);
}

static void contour_insertion(BMesh *bm, BMesh *bm_orig, BLI_Buffer *new_vert_buffer,
							  BLI_Buffer *cusp_edges, struct OpenSubdiv_EvaluatorDescr *eval,
							  const float cam_loc[3]){
    int i, face_index;
	BMEdge *e;
	BMFace *f;
	BMIter iter_e;
	float v1_u, v1_v, v2_u, v2_v;

	int orig_verts = BM_mesh_elem_count(bm_orig, BM_VERT);
	int orig_edges = BM_mesh_elem_count(bm_orig, BM_EDGE);
    int initial_edges = BM_mesh_elem_count(bm, BM_EDGE);

	BM_ITER_MESH_INDEX (e, &iter_e, bm, BM_EDGES_OF_MESH, i) {
        if( !(i < initial_edges) ){
			//Now we are working on edges we added in this function
			break;
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
				v_buf1 = BLI_buffer_at(new_vert_buffer, Vert_buf, v1_idx + 1 - orig_verts);
				v1_u = v_buf1.u;
				v1_v = v_buf2.v;
				if( v_buf1.orig_face ){
					v1_has_face = true;
				}
			} else {
				v1 = BM_vert_at_index_find( bm_orig, BM_elem_index_get( e->v1 ) ); 
			}
			if( (v2_idx + 1) > orig_verts){
				v_buf2 = BLI_buffer_at(new_vert_buffer, Vert_buf, v2_idx + 1 - orig_verts);
				v2_u = v_buf2.u;
				v2_v = v_buf2.v;
				if( v_buf2.orig_face ){
					v2_has_face = true;
				}
			} else {
				v2 = BM_vert_at_index_find( bm_orig, BM_elem_index_get( e->v2 ) ); 
			}

			if( v1 && v2 ){
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
				if( v2_has_face ){
					f = v_buf2.orig_face; 
				} else {
					BMVert *vert_arr[3];

                    vert_arr[0] = v1;
					vert_arr[1] = v_buf2.orig_edge->v1;
					vert_arr[2] = v_buf2.orig_edge->v2;
                    //TODO check if this fails
					BM_face_exists_overlap(vert_arr, 3, &f);
					convert_uv_to_new_face( v_buf2.orig_edge, f, &v2_u, &v2_v);
				}
				get_uv_coord(v1, f, &v1_u, &v1_v);
			} else if ( v2 ){
				if( v1_has_face ){
					f = v_buf1.orig_face; 
				} else {
					BMVert *vert_arr[3];

                    vert_arr[0] = v2;
					vert_arr[1] = v_buf1.orig_edge->v1;
					vert_arr[2] = v_buf1.orig_edge->v2;
                    //TODO check if this fails
					BM_face_exists_overlap(vert_arr, 3, &f);
					convert_uv_to_new_face( v_buf1.orig_edge, f, &v1_u, &v1_v);
				}
				get_uv_coord(v2, f, &v2_u, &v2_v);
			} else {
				//TODO This should not happen. Check if this really is the case.
                printf("Two verts on the same orig face\n");
				f = v_buf1.orig_face;
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

static struct OpenSubdiv_EvaluatorDescr *create_osd_eval(BMesh *bm){
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

    //subsurf_make_derived_from_derived( 
	
	//opensubdiv_ensureEvaluator(  <<--- USE this
	
	//opensubdiv_initEvaluatorFace(
	
	//opensubdiv_initEvaluator(

	//ccgSubSurf_free( <<--- for freeing "ss"

	DerivedMesh *result;
	BMesh *bm_orig, *bm;

	//CCGSubSurf *ss;
    struct OpenSubdiv_EvaluatorDescr *osd_eval;
    
	//ss = get_ss_for_osd(dm);
	
	//TODO fix the get_osd_eval declaration ("WITH_OPENSUBDIV" problem)
	//osd_eval = get_osd_eval(ss);

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

		split_BB_FF_edges(bm, bm_orig, osd_eval, cam_loc, &new_vert_buffer);

		// (6.2) Contour Insertion

		contour_insertion(bm, bm_orig, &new_vert_buffer, NULL, osd_eval, cam_loc);

		debug_colorize(bm, cam_loc);
		BLI_buffer_free(&new_vert_buffer);
	}
	result = CDDM_from_bmesh(bm, true);

	//ccgSubSurf_free(ss);
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
