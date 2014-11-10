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

void verts_to_limit(BMesh *bm, CCGSubSurf *ss){
	
	int i;
	float new_co[3];
	
	BMIter iter_v, iter_f;
	BMVert *v;
    BMFace *f;
	BM_ITER_MESH_INDEX (f, &iter_f, bm, BM_FACES_OF_MESH, i) {
			BM_ITER_ELEM (v, &iter_v, f, BM_VERTS_OF_FACE) {
				new_co = {0.0f, 0.0f, 0.0f};
				openSubdiv_evaluateLimit(ss->osd_evaluator, i, 0, 0, new_co, NULL, NULL);
				v.co = new_co;
			}
	}

}

/* bmesh only function */
static DerivedMesh *mybmesh_do(DerivedMesh *dm, MyBMeshModifierData *mmd)
{

    //subsurf_make_derived_from_derived( 
	
	//opensubdiv_ensureEvaluator(  <<--- USE this
	
	//opensubdiv_initEvaluatorFace(
	
	//opensubdiv_initEvaluator(

	//ccgSubSurf_free( <<--- for freeing "ss"

    //From subsurf_calculate_limit_positions
	CCGSubSurf *ss = _getSubSurf(NULL, 1, 3, CCG_USE_ARENA);
	ss_sync_from_derivedmesh(ss, dm, NULL, 0);
	
    //Get opensubdiv evaluator inited
	opensubdiv_ensureEvaluator(ss);
	opensubdiv_updateCoarsePositions(ss);

	DerivedMesh *result;
	BMesh *bm_orig, *bm;

	bm = DM_to_bmesh(dm, true);
    
    verts_to_limit(bm, ss);

	BM_mesh_triangulate(bm, MOD_TRIANGULATE_QUAD_FIXED, MOD_TRIANGULATE_NGON_BEAUTY, false, NULL, NULL);

	result = CDDM_from_bmesh(bm, true);


	ccgSubSurf_free(ss);
	BM_mesh_free(bm);

	result->dirty |= DM_DIRTY_NORMALS;

	//TODO use this to check if we need to subdivide the mesh to get a quad mesh.
	if(0){
	BMIter iter;
	BMFace *f, *f_next;
	/* use the mutable iterator so we can remove data as its looped over */
	BM_ITER_MESH_MUTYABLE (f, f_next, &iter, bm, BM_FACES_OF_MESH) {
		if (f->len == mmd->sides) {
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

	mmd->sides = 4;
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

static DerivedMesh *applyModifier(ModifierData *md, struct Object *UNUSED(ob),
                                  DerivedMesh *dm,
                                  ModifierApplyFlag UNUSED(flag))
{
	DerivedMesh *result;

	MyBMeshModifierData *mmd = (MyBMeshModifierData *)md;

	if (!(result = triangulate_dm(dm, mmd))) {
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
	/* updateDepgraph */    NULL,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */  dependsOnNormals,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};
