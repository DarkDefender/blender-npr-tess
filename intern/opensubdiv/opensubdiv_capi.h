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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2013 Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Sergey Sharybin
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __OPENSUBDIV_CAPI_H__
#define __OPENSUBDIV_CAPI_H__

#undef OPENSUBDIV_LEGACY_DRAW

#ifdef __cplusplus
extern "C" {
#endif

// Types declaration.
struct OpenSubdiv_EvaluatorDescr;
struct OpenSubdiv_GLMesh;

typedef struct OpenSubdiv_EvaluatorDescr OpenSubdiv_EvaluatorDescr;
typedef struct OpenSubdiv_GLMesh OpenSubdiv_GLMesh;

#ifdef __cplusplus
struct OpenSubdiv_GLMeshDescr;
typedef struct OpenSubdiv_GLMesh {
	int controller_type;
	OpenSubdiv_GLMeshDescr *descriptor;
	int level;
} OpenSubdiv_GLMesh;
#endif

// Keep this a bitmask os it's possible to pass vailable
// controllers to Blender.
enum {
	OPENSUBDIV_CONTROLLER_CPU                      = (1 << 0),
	OPENSUBDIV_CONTROLLER_OPENMP                   = (1 << 1),
	OPENSUBDIV_CONTROLLER_OPENCL                   = (1 << 2),
	OPENSUBDIV_CONTROLLER_CUDA                     = (1 << 3),
	OPENSUBDIV_CONTROLLER_GLSL_TRANSFORM_FEEDBACK  = (1 << 4),
	OPENSUBDIV_CONTROLLER_GLSL_COMPUTE             = (1 << 5),
};

enum {
	OPENSUBDIV_SCHEME_CATMARK,
	OPENSUBDIV_SCHEME_BILINEAR,
	OPENSUBDIV_SCHEME_LOOP,
};

OpenSubdiv_GLMesh *openSubdiv_createOsdGLMeshFromEvaluator(
    OpenSubdiv_EvaluatorDescr *evaluator_descr,
    int controller_type,
    int level,
    int scheme);

void openSubdiv_deleteOsdGLMesh(OpenSubdiv_GLMesh *gl_mesh);
unsigned int openSubdiv_getOsdGLMeshPatchIndexBuffer(
        OpenSubdiv_GLMesh *gl_mesh);
unsigned int openSubdiv_getOsdGLMeshVertexBuffer(OpenSubdiv_GLMesh *gl_mesh);
void openSubdiv_osdGLMeshUpdateVertexBuffer(OpenSubdiv_GLMesh *gl_mesh,
                                            const float *vertex_data,
                                            int start_vertex,
                                            int num_verts);
void openSubdiv_osdGLMeshRefine(OpenSubdiv_GLMesh *gl_mesh);
void openSubdiv_osdGLMeshSynchronize(OpenSubdiv_GLMesh *gl_mesh);
void openSubdiv_osdGLMeshBindVertexBuffer(OpenSubdiv_GLMesh *gl_mesh);

/* ** Initialize/Deinitialize global OpenGL drawing buffers/GLSL programs ** */
void openSubdiv_osdGLDisplayInit(void);
void openSubdiv_osdGLDisplayDeinit(void);

/* ** Actual drawing ** */

/* Initialize all the invariants which stays the same for every single path,
 * for example lighting model stays untouched for the whole mesh.
 *
 * TODO(sergey): Some of the stuff could be initialized once for all meshes.
 */
void openSubdiv_osdGLMeshDisplayPrepare(int use_osd_glsl);

/* Draw patches which corresponds to a given partition. */
void openSubdiv_osdGLMeshDisplay(OpenSubdiv_GLMesh *gl_mesh,
                                 int fill_quads,
                                 int start_partition,
                                 int num_partitions);

/* ** Utility functions ** */
int openSubdiv_supportGPUDisplay(void);
int openSubdiv_getAvailableControllers(void);
void openSubdiv_cleanup(void);

#ifdef __cplusplus
}
#endif

#endif  // __OPENSUBDIV_CAPI_H__
