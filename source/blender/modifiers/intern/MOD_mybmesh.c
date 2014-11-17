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

//TODO for Kr look in subdiv.cpp in coutours source code (II)

static void verts_to_limit(BMesh *bm, struct OpenSubdiv_EvaluatorDescr *eval){
	
	int i, j;
	
	BMIter iter_v, iter_f;
	BMVert *vert;
    BMFace *f;

    //TODO is it possible to only get non adjusted verts?
	//IE not moving a vert more than once.

	BM_ITER_MESH_INDEX (f, &iter_f, bm, BM_FACES_OF_MESH, i) {
			BM_ITER_ELEM_INDEX (vert, &iter_v, f, BM_VERTS_OF_FACE, j) {
				float u,v;
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
				openSubdiv_evaluateLimit(eval, i, u, v, vert->co, NULL, NULL);
				//printf("j: %d\n",j);
			}
			//printf("i: %d\n",i);
			//printf("face i: %d\n", BM_elem_index_get(f));
	}

}

static bool calc_if_F(const float cam_loc[3], const float P[3], const float du[3], const float dv[3]){
	//Is the point front facing?
	float nor[3], view_vec[3];

	cross_v3_v3v3(nor, dv, du);
	//TODO normalization is probably not needed
	normalize_v3(nor);
	sub_v3_v3v3(view_vec, cam_loc, P);

	return ( dot_v3v3(nor, view_vec) < 0);
}

static void split_edge_and_move_vert(BMesh *bm, BMEdge *edge, float new_pos[3]){
	//Split edge one time and move the created vert to new_pos
	
    BMVert *v;
	printf("Split edge!\n");

	BM_mesh_elem_hflag_disable_all(bm, BM_EDGE, BM_ELEM_TAG, false);
	BM_elem_flag_enable(edge, BM_ELEM_TAG);
	BMO_op_callf(bm, BMO_FLAG_DEFAULTS,
			"subdivide_edges edges=%he cuts=%i quad_corner_type=%i use_single_edge=true",
			BM_ELEM_TAG, 1, SUBD_STRAIGHT_CUT);
	//Get the newly created vertex
	v = BM_vert_at_index_find(bm, BM_mesh_elem_count(bm, BM_VERT) - 1);
	copy_v3_v3(v->co, new_pos);
	//TODO save what face the new vert belonged to. Also u,v coords
}

static void split_BB_FF_edges(BMesh *bm, BMesh *bm_orig, struct OpenSubdiv_EvaluatorDescr *eval, const float cam_loc[3]){
    //Split BB,FF edges if they have sign crossings
	
	int i, j, face_index;
	BMIter iter_f, iter_e, iter_v;
	BMEdge *e;
	BMFace *f;
	BMVert *vert, *v1, *v2;
	float v1_u, v1_v, v2_u, v2_v, step;
	bool is_F;
	int orig_edges = BM_mesh_elem_count(bm_orig, BM_EDGE);
    int initial_edges = BM_mesh_elem_count(bm, BM_EDGE);

	//Do 10 samples but don't check end and start point
	step = 1.0f/12.0f;

	BM_ITER_MESH_INDEX (e, &iter_e, bm, BM_EDGES_OF_MESH, i) {
        if( !(i < initial_edges) ){
			//Now we are working on edges we added in this function
			break;
		}

		if( i < orig_edges ){
            //This edge exists on the original mesh
			//TODO why do I have to use find? Segfault otherwise...
			//remember to replace the rest of "at_index"
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
		} else {
			BMFace *f2;
			BMIter iter_f2;
			bool found_face = false;

			//This should be safe because the vert count is still the same as the original mesh.
			v1 = BM_vert_at_index_find( bm_orig, BM_elem_index_get( e->v1 ) ); 
			v2 = BM_vert_at_index_find( bm_orig, BM_elem_index_get( e->v2 ) ); 

			//TODO is there a better way to find the shared face?
			BM_ITER_ELEM (f, &iter_f, v1, BM_FACES_OF_VERT) {
				BM_ITER_ELEM (f2, &iter_f2, v2, BM_FACES_OF_VERT) {
					if( f == f2 ){
						found_face = true;
						face_index = BM_elem_index_get(f);
						break;
					}
				}
				if(found_face){
					break;
				}
			}
		}

		BM_ITER_ELEM_INDEX (vert, &iter_v, f, BM_VERTS_OF_FACE, j) {
			if(v1 == vert){
				switch(j){
					case 1 : v1_u = 1, v1_v = 0;
							 break;
					case 2 : v1_u = v1_v = 1;
							 break;
					case 3 : v1_u = 0, v1_v = 1;
							 break;
					default: v1_u = v1_v = 0;
							 break;
				}
			} else if(v2 == vert){
				switch(j){
					case 1 : v2_u = 1, v2_v = 0;
							 break;
					case 2 : v2_u = v2_v = 1;
							 break;
					case 3 : v2_u = 0, v2_v = 1;
							 break;
					default: v2_u = v2_v = 0;
							 break;
				}
			}	
				
		}

		{
			float P1[3], P2[3], du[3], dv[3];
			//is this a FF or BB edge?
			openSubdiv_evaluateLimit(eval, face_index, v1_u, v1_v, P1, du, dv);

			is_F = calc_if_F(cam_loc, P1, du, dv);

			openSubdiv_evaluateLimit(eval, face_index, v2_u, v2_v, P2, du, dv);

			if( is_F  != calc_if_F(cam_loc, P2, du, dv) ){
				//FB edge, we only want to split FF or BB
				//Skip to next edge
				continue;  
			}
		}

		{
            int i;
            float u, v;
			float P[3], du[3], dv[3];

			if( v1_u == v2_u ){
				u = v1_u;
				v = step;

				for(i=0; i < 10; i++){
					v += step;
					openSubdiv_evaluateLimit(eval, face_index, u, v, P, du, dv);
					if( calc_if_F(cam_loc, P, du, dv) != is_F ){
						split_edge_and_move_vert(bm, e, P);
						break;
					}
				}
			} else if ( v1_v == v2_v ){
				u = step;
				v = v1_v;

				for(i=0; i < 10; i++){
					u += step;
					openSubdiv_evaluateLimit(eval, face_index, u, v, P, du, dv);
					if( calc_if_F(cam_loc, P, du, dv) != is_F ){
						split_edge_and_move_vert(bm, e, P);
						break;
					}
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
					u += step_u;
					v += step;
					openSubdiv_evaluateLimit(eval, face_index, u, v, P, du, dv);
					if( calc_if_F(cam_loc, P, du, dv) != is_F ){
						split_edge_and_move_vert(bm, e, P);
						break;
					}
				}
			}
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

	split_BB_FF_edges(bm, bm_orig, osd_eval, cam_loc);

	// (6.2) Contour Insertion

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
		float tmp[4][4];

		invert_m4_m4(tmp, mmd->camera_ob->obmat);
		copy_v3_v3(cam_loc, mmd->camera_ob->loc);
		//convert camera origin to world coord and the to local modifier obj coords
		mul_m4_v3(tmp, cam_loc);
		mul_m4_v3(ob->obmat, cam_loc);
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
