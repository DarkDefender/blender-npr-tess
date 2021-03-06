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
 * Contributor(s): Sergey Sharybin.
 *                 Brecht van Lommel
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#include "opensubdiv_capi.h"

#ifdef _MSC_VER
#  include "iso646.h"
#endif

#include <GL/glew.h>

#include <opensubdiv/osd/glMesh.h>
#include <opensubdiv/osd/glDrawRegistry.h>
#include <opensubdiv/osdutil/evaluator_capi.h>
#include <opensubdiv/osdutil/topology.h>
#include <opensubdiv/osdutil/mesh.h>

// CPU Backend
#include <opensubdiv/osd/cpuGLVertexBuffer.h>
#include <opensubdiv/osd/cpuComputeContext.h>
#include <opensubdiv/osd/cpuComputeController.h>

#ifdef OPENSUBDIV_HAS_OPENMP
#  include <opensubdiv/osd/ompComputeController.h>
#endif

#ifdef OPENSUBDIV_HAS_OPENCL
#  include <opensubdiv/osd/clGLVertexBuffer.h>
#  include <opensubdiv/osd/clComputeContext.h>
#  include <opensubdiv/osd/clComputeController.h>
#  include "clInit.h"
#endif

#ifdef OPENSUBDIV_HAS_CUDA
#  include <opensubdiv/osd/cudaGLVertexBuffer.h>
#  include <opensubdiv/osd/cudaComputeContext.h>
#  include <opensubdiv/osd/cudaComputeController.h>
#  include "cudaInit.h"
#endif

#ifdef OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK
#  include <opensubdiv/osd/glslTransformFeedbackComputeContext.h>
#  include <opensubdiv/osd/glslTransformFeedbackComputeController.h>
#  include <opensubdiv/osd/glVertexBuffer.h>
#endif

#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
#  include <opensubdiv/osd/glslComputeContext.h>
#  include <opensubdiv/osd/glslComputeController.h>
#  include <opensubdiv/osd/glVertexBuffer.h>
#endif

#include <opensubdiv/osdutil/patchPartitioner.h>

#include "opensubdiv_partitioned.h"

#include "MEM_guardedalloc.h"

// **************** Types declaration ****************

struct OpenSubdiv_ComputeControllerDescr;

typedef struct OpenSubdiv_ComputeController {
	int type;
	OpenSubdiv_ComputeControllerDescr *descriptor;
} OpenSubdiv_ComputeController;

using OpenSubdiv::OsdCpuComputeController;
using OpenSubdiv::OsdGLDrawContext;
using OpenSubdiv::OsdGLMeshInterface;
using OpenSubdiv::OsdMesh;
using OpenSubdiv::OsdMeshBitset;
using OpenSubdiv::OsdUtilSubdivTopology;
using OpenSubdiv::OsdVertex;

typedef OpenSubdiv::HbrMesh<OsdVertex> OsdHbrMesh;
typedef OpenSubdiv::HbrFace<OpenSubdiv::OsdVertex> OsdHbrFace;

#if defined(OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK) || \
    defined(OPENSUBDIV_HAS_GLSL_COMPUTE)
using OpenSubdiv::OsdGLVertexBuffer;
#endif

using OpenSubdiv::OsdGLDrawContext;
using OpenSubdiv::PartitionedMesh;

// CPU backend
using OpenSubdiv::OsdCpuGLVertexBuffer;
using OpenSubdiv::OsdCpuComputeController;
static OpenSubdiv_ComputeController *g_cpuComputeController = NULL;
typedef PartitionedMesh<OsdCpuGLVertexBuffer,
                        OsdCpuComputeController,
                        OsdGLDrawContext> OsdCpuMesh;

#ifdef OPENSUBDIV_HAS_OPENMP
using OpenSubdiv::OsdOmpComputeController;
static OpenSubdiv_ComputeController *g_ompComputeController = NULL;
typedef PartitionedMesh<OsdCpuGLVertexBuffer,
                        OsdOmpComputeController,
                        OsdGLDrawContext> OsdOmpMesh;
#endif

#ifdef OPENSUBDIV_HAS_OPENCL
using OpenSubdiv::OsdCLGLVertexBuffer;
using OpenSubdiv::OsdCLComputeController;
typedef PartitionedMesh<OsdCLGLVertexBuffer,
                        OsdCLComputeController,
                        OsdGLDrawContext> OsdCLMesh;
static OpenSubdiv_ComputeController *g_clComputeController = NULL;
static cl_context g_clContext;
static cl_command_queue g_clQueue;
#endif

#ifdef OPENSUBDIV_HAS_CUDA
using OpenSubdiv::OsdCudaComputeController;
using OpenSubdiv::OsdCudaGLVertexBuffer;
typedef PartitionedMesh<OsdCudaGLVertexBuffer,
                        OsdCudaComputeController,
                        OsdGLDrawContext> OsdCudaMesh;
static OpenSubdiv_ComputeController *g_cudaComputeController = NULL;
static bool g_cudaInitialized = false;
#endif

#ifdef OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK
using OpenSubdiv::OsdGLSLTransformFeedbackComputeController;
typedef PartitionedMesh<OsdGLVertexBuffer,
                        OsdGLSLTransformFeedbackComputeController,
                        OsdGLDrawContext> OsdGLSLTransformFeedbackMesh;
static OpenSubdiv_ComputeController *g_glslTransformFeedbackComputeController = NULL;
#endif

#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
using OpenSubdiv::OsdGLSLComputeController;
typedef PartitionedMesh<OsdGLVertexBuffer,
                        OsdGLSLComputeController,
                        OsdGLDrawContext> OsdGLSLComputeMesh;
static OpenSubdiv_ComputeController *g_glslOComputeComputeController = NULL;
#endif

static OpenSubdiv_ComputeController *alloc_controller(int controller_type)
{
	OpenSubdiv_ComputeController *controller =
		(OpenSubdiv_ComputeController *) OBJECT_GUARDED_NEW(
			OpenSubdiv_ComputeController);
	controller->type = controller_type;
	return controller;
}

static OpenSubdiv_ComputeController *openSubdiv_getController(
	int controller_type)
{
#ifdef OPENSUBDIV_HAS_OPENCL
	if (controller_type == OPENSUBDIV_CONTROLLER_OPENCL &&
	    g_clContext == NULL)
	{
		if (!HAS_CL_VERSION_1_1()) {
			printf("OpenCL is not supported on this system\n");
			return NULL;
		}
		if (initCL(&g_clContext, &g_clQueue) == false) {
			printf("Error in initializing OpenCL\n");
			return NULL;
		}
	}
#endif

#ifdef OPENSUBDIV_HAS_CUDA
	if (controller_type == OPENSUBDIV_CONTROLLER_CUDA &&
	   g_cudaInitialized == false)
	{
		if (!HAS_CUDA_VERSION_4_0()) {
			printf("CUDA is not supported on this system\n");
			return NULL;
		}
		g_cudaInitialized = true;
	}
#endif

	switch (controller_type) {
#define CHECK_CONTROLLER_TYPR(type, var, class) \
	case OPENSUBDIV_CONTROLLER_ ## type: \
		if (var == NULL) { \
			var = alloc_controller(controller_type); \
			var->descriptor = \
				(OpenSubdiv_ComputeControllerDescr *) new class(); \
		} \
		return var;

		CHECK_CONTROLLER_TYPR(CPU,
		                      g_cpuComputeController,
		                      OsdCpuComputeController);

#ifdef OPENSUBDIV_HAS_OPENMP
		CHECK_CONTROLLER_TYPR(OPENMP,
		                      g_ompComputeController,
		                      OsdOmpComputeController);
#endif

#ifdef OPENSUBDIV_HAS_OPENCL
	case OPENSUBDIV_CONTROLLER_OPENCL:
		if (g_clComputeController == NULL) {
			g_clComputeController = alloc_controller(controller_type);
			g_clComputeController->descriptor =
				(OpenSubdiv_ComputeControllerDescr *)
					new OsdCLComputeController(g_clContext, g_clQueue);
		}
		return g_clComputeController;
#endif

#ifdef OPENSUBDIV_HAS_CUDA
		CHECK_CONTROLLER_TYPR(CUDA,
		                      g_cudaComputeController,
		                      OsdCudaComputeController);
#endif

#ifdef OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK
		CHECK_CONTROLLER_TYPR(GLSL_TRANSFORM_FEEDBACK,
		                      g_glslTransformFeedbackComputeController,
		                      OsdGLSLTransformFeedbackComputeController);
#endif

#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
		CHECK_CONTROLLER_TYPR(GLSL_COMPUTE,
		                      g_glslComputeController,
		                      OsdGLSLComputeController);
#endif
	}

	return NULL;
}

static OpenSubdiv::OsdUtilMesh<OsdVertex>::Scheme get_osd_scheme(int scheme)
{
	switch (scheme) {
		case OPENSUBDIV_SCHEME_CATMARK:
			return OpenSubdiv::OsdUtilMesh<OsdVertex>::SCHEME_CATMARK;
		case OPENSUBDIV_SCHEME_BILINEAR:
			return OpenSubdiv::OsdUtilMesh<OsdVertex>::SCHEME_BILINEAR;
		case OPENSUBDIV_SCHEME_LOOP:
			return OpenSubdiv::OsdUtilMesh<OsdVertex>::SCHEME_LOOP;
		default:
			assert(!"Wrong subdivision scheme");
	}
	return OpenSubdiv::OsdUtilMesh<OsdVertex>::SCHEME_BILINEAR;
}

/* TODO(sergey): Currently we use single coarse face per partition,
 * which allows to have per-face material assignment but which also
 * increases number of glDrawElements() calls.
 *
 * Ideally here we need to partition like this, but do some conjunction
 * at draw time, so adjacent faces with the same material are displayed
 * in a single chunk.
 */
static void get_partition_per_face(OsdHbrMesh &hmesh,
                                   std::vector<int> *idsOnPtexFaces)
{
	int numFaces = hmesh.GetNumCoarseFaces();

	// First, assign partition ID to each coarse face.
	std::vector<int> idsOnCoarseFaces;
	for (int i = 0; i < numFaces; ++i) {
		int partitionID = i;
		idsOnCoarseFaces.push_back(partitionID);
	}

	// Create ptex index to coarse face index mapping.
	OsdHbrFace *lastFace = hmesh.GetFace(numFaces - 1);
	int numPtexFaces = lastFace->GetPtexIndex();
	numPtexFaces += (hmesh.GetSubdivision()->FaceIsExtraordinary(&hmesh,
	                                                             lastFace) ?
	                 lastFace->GetNumVertices() : 1);

	// TODO(sergey): Duplicated logic to simpleHbr.
	std::vector<int> ptexIndexToFaceMapping(numPtexFaces);
	int ptexIndex = 0;
	for (int i = 0; i < numFaces; ++i) {
		OsdHbrFace *f = hmesh.GetFace(i);
		ptexIndexToFaceMapping[ptexIndex++] = i;
		int numVerts = f->GetNumVertices();
		if (numVerts != 4 ) {
			for (int j = 0; j < numVerts-1; ++j) {
				ptexIndexToFaceMapping[ptexIndex++] = i;
			}
		}
	}
	assert((int)ptexIndexToFaceMapping.size() == numPtexFaces);

	// Convert ID array from coarse face index space to ptex index space.
	for (int i = 0; i < numPtexFaces; ++i) {
		idsOnPtexFaces->push_back(idsOnCoarseFaces[ptexIndexToFaceMapping[i]]);
	}
}

/* TODO(sergey): It might be good to have in evaluator_capi. */
void openSubdiv_evlauatorClearTags(
    OpenSubdiv_EvaluatorDescr *evaluator_descr)
{
	OsdUtilSubdivTopology *topology =
	    (OsdUtilSubdivTopology *)openSubdiv_getEvaluatorTopologyDescr(
	        evaluator_descr);

	topology->tagData.tags.clear();
	topology->tagData.intArgs.clear();
	topology->tagData.floatArgs.clear();
	topology->tagData.stringArgs.clear();
}

void openSubdiv_evaluatorSetEdgeSharpness(
    OpenSubdiv_EvaluatorDescr *evaluator_descr,
    int v0, int v1,
    float sharpness)
{
	OsdUtilSubdivTopology *topology =
	    (OsdUtilSubdivTopology *)openSubdiv_getEvaluatorTopologyDescr(
	        evaluator_descr);
	int indices[] = {v0, v1};

	topology->tagData.AddCrease(indices, 2, &sharpness, 1);
}

const float *openSubdiv_evaluatorGetFloatTagArgs(
    OpenSubdiv_EvaluatorDescr *evaluator_descr)
{
	OsdUtilSubdivTopology *topology =
	    (OsdUtilSubdivTopology *)openSubdiv_getEvaluatorTopologyDescr(
	        evaluator_descr);
	return &topology->tagData.floatArgs[0];
}

struct OpenSubdiv_GLMesh *openSubdiv_createOsdGLMeshFromEvaluator(
    OpenSubdiv_EvaluatorDescr *evaluator_descr,
    int controller_type,
    int level,
    int scheme,
    int subdivide_uvs)
{
	OpenSubdiv_ComputeController *controller =
		openSubdiv_getController(controller_type);
	if (controller == NULL) {
		return NULL;
	}

	OsdUtilSubdivTopology *topology;
	OpenSubdiv::OsdUtilMesh<OsdVertex> util_mesh;

	topology = (OsdUtilSubdivTopology *)openSubdiv_getEvaluatorTopologyDescr(
		evaluator_descr);

	if (util_mesh.Initialize(*topology,
	                         NULL,
	                         get_osd_scheme(scheme)) == false) {
		return NULL;
	}

	OsdHbrMesh *hmesh = util_mesh.GetHbrMesh();

	std::vector<int> idsOnPtexFaces;
	get_partition_per_face(*hmesh, &idsOnPtexFaces);

	hmesh->SetFVarInterpolateBoundaryMethod(
	        subdivide_uvs
	            ? OsdHbrMesh::k_InterpolateBoundaryEdgeAndCorner
	            : OsdHbrMesh::k_InterpolateBoundaryNone);
	hmesh->SetFVarPropagateCorners(false);

	OsdMeshBitset bits;
	/* TODO(sergey): Adaptive subdivisions are not currently
	 * possible because of the lack of tessellation shader.
	 */
	bits.set(OpenSubdiv::MeshAdaptive, 0);
	bits.set(OpenSubdiv::MeshFVarData, 1);

	int num_vertex_elements = 6;
	int num_varying_elements = 0;

	OsdGLMeshInterface *mesh = NULL;

	switch (controller_type) {
#define CHECK_CONTROLLER_TYPE(type, class, controller_class) \
		case OPENSUBDIV_CONTROLLER_ ## type: \
			mesh = (OsdGLMeshInterface *) \
				new class( \
					(controller_class *) controller->descriptor, \
					hmesh, \
					num_vertex_elements, \
					num_varying_elements, \
					level, \
					bits, \
					idsOnPtexFaces); \
			break;
		CHECK_CONTROLLER_TYPE(CPU, OsdCpuMesh, OsdCpuComputeController)

#ifdef OPENSUBDIV_HAS_OPENMP
		CHECK_CONTROLLER_TYPE(OPENMP, OsdOmpMesh, OsdOmpComputeController)
#endif

#ifdef OPENSUBDIV_HAS_OPENCL
		case OPENSUBDIV_CONTROLLER_OPENCL:
			mesh = (OsdGLMeshInterface *)
				new OsdCLMesh(
					(OsdCLComputeController *) controller->descriptor,
					hmesh,
					num_vertex_elements,
					num_varying_elements,
					level,
					bits,
					g_clContext,
					g_clQueue,
					idsOnPtexFaces);
			break;
#endif

#ifdef OPENSUBDIV_HAS_CUDA
		CHECK_CONTROLLER_TYPE(CUDA, OsdCudaMesh, OsdCudaComputeController)
#endif

#ifdef OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK
		CHECK_CONTROLLER_TYPE(GLSL_TRANSFORM_FEEDBACK,
		                      OsdGLSLTransformFeedbackMesh,
		                      OsdGLSLTransformFeedbackComputeController)
#endif

#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
		CHECK_CONTROLLER_TYPE(GLSL_COMPUTE,
		                      OsdGLSLComputeMesh,
		                      OsdGLSLComputeController)
#endif

#undef CHECK_CONTROLLER_TYPE
	}

	if (mesh == NULL) {
		return NULL;
	}

	OpenSubdiv_GLMesh *gl_mesh =
		(OpenSubdiv_GLMesh *) OBJECT_GUARDED_NEW(OpenSubdiv_GLMesh);
	gl_mesh->controller_type = controller_type;
	gl_mesh->descriptor = (OpenSubdiv_GLMeshDescr *) mesh;
	gl_mesh->level = level;

	return gl_mesh;
}

void openSubdiv_deleteOsdGLMesh(struct OpenSubdiv_GLMesh *gl_mesh)
{
	switch (gl_mesh->controller_type) {
#define CHECK_CONTROLLER_TYPE(type, class) \
		case OPENSUBDIV_CONTROLLER_ ## type: \
			delete (class *) gl_mesh->descriptor; \
			break;

		CHECK_CONTROLLER_TYPE(CPU, OsdCpuMesh)

#ifdef OPENSUBDIV_HAS_OPENMP
		CHECK_CONTROLLER_TYPE(OPENMP, OsdOmpMesh)
#endif

#ifdef OPENSUBDIV_HAS_OPENCL
		CHECK_CONTROLLER_TYPE(OPENCL, OsdCLMesh)
#endif

#ifdef OPENSUBDIV_HAS_CUDA
		CHECK_CONTROLLER_TYPE(CUDA, OsdCudaMesh)
#endif

#ifdef OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK
		CHECK_CONTROLLER_TYPE(GLSL_TRANSFORM_FEEDBACK,
		                      OsdGLSLTransformFeedbackMesh)
#endif

#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
		CHECK_CONTROLLER_TYPE(GLSL_COMPUTE, OsdGLSLComputeMesh)
#endif

#undef CHECK_CONTROLLER_TYPE
	}

	OBJECT_GUARDED_DELETE(gl_mesh, OpenSubdiv_GLMesh);
}

unsigned int openSubdiv_getOsdGLMeshPatchIndexBuffer(struct OpenSubdiv_GLMesh *gl_mesh)
{
	return ((OsdGLMeshInterface *)gl_mesh->descriptor)->GetDrawContext()->GetPatchIndexBuffer();
}

unsigned int openSubdiv_getOsdGLMeshVertexBuffer(struct OpenSubdiv_GLMesh *gl_mesh)
{
	return ((OsdGLMeshInterface *)gl_mesh->descriptor)->BindVertexBuffer();
}

void openSubdiv_osdGLMeshUpdateVertexBuffer(struct OpenSubdiv_GLMesh *gl_mesh,
                                            const float *vertex_data,
                                            int start_vertex,
                                            int num_verts)
{
	((OsdGLMeshInterface *)gl_mesh->descriptor)->UpdateVertexBuffer(vertex_data,
	                                                                start_vertex,
	                                                                num_verts);
}

void openSubdiv_osdGLMeshRefine(struct OpenSubdiv_GLMesh *gl_mesh)
{
	((OsdGLMeshInterface *)gl_mesh->descriptor)->Refine();
}

void openSubdiv_osdGLMeshSynchronize(struct OpenSubdiv_GLMesh *gl_mesh)
{
	((OsdGLMeshInterface *)gl_mesh->descriptor)->Synchronize();
}

void openSubdiv_osdGLMeshBindVertexBuffer(OpenSubdiv_GLMesh *gl_mesh)
{
	((OsdGLMeshInterface *)gl_mesh->descriptor)->BindVertexBuffer();
}

int openSubdiv_supportGPUDisplay(void)
{
	return GL_EXT_geometry_shader4 &&
	       GL_ARB_gpu_shader5 &&
	       glProgramParameteriEXT;
}

int openSubdiv_getAvailableControllers(void)
{
	if (!openSubdiv_supportGPUDisplay()) {
		return 0;
	}

	int flags = OPENSUBDIV_CONTROLLER_CPU;

#ifdef OPENSUBDIV_HAS_OPENMP
	flags |= OPENSUBDIV_CONTROLLER_OPENMP;
#endif

#ifdef OPENSUBDIV_HAS_OPENCL
	if (HAS_CL_VERSION_1_1()) {
		flags |= OPENSUBDIV_CONTROLLER_OPENCL;
	}
#endif

#ifdef OPENSUBDIV_HAS_CUDA
	if (HAS_CUDA_VERSION_4_0()) {
		flags |= OPENSUBDIV_CONTROLLER_CUDA;
	}
#endif

#ifdef OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK
	if (GLEW_ARB_texture_buffer_object) {
		flags |= OPENSUBDIV_CONTROLLER_GLSL_TRANSFORM_FEEDBACK;
	}
#endif

#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
	flags |= OPENSUBDIV_CONTROLLER_GLSL_COMPUTE;
#endif

	return flags;
}

void openSubdiv_cleanup(void)
{
	openSubdiv_osdGLDisplayDeinit();

#define DELETE_DESCRIPTOR(var, class) \
	if (var != NULL) { \
		delete (class*) var->descriptor; \
		OBJECT_GUARDED_DELETE(var, OpenSubdiv_ComputeController); \
	}

	DELETE_DESCRIPTOR(g_cpuComputeController,
	                  OsdCpuComputeController);

#ifdef OPENSUBDIV_HAS_OPENMP
	DELETE_DESCRIPTOR(g_ompComputeController,
	                  OsdOmpComputeController);
#endif

#ifdef OPENSUBDIV_HAS_OPENCL
	DELETE_DESCRIPTOR(g_clComputeController,
	                  OsdCLComputeController);
    uninitCL(g_clContext, g_clQueue);
#endif

#ifdef OPENSUBDIV_HAS_CUDA
	DELETE_DESCRIPTOR(g_cudaComputeController,
	                  OsdCudaComputeController);
#endif

#ifdef OPENSUBDIV_HAS_GLSL_TRANSFORM_FEEDBACK
	if (GLEW_ARB_texture_buffer_object) {
		DELETE_DESCRIPTOR(g_glslTransformFeedbackComputeController,
		                  OsdGLSLTransformFeedbackComputeController);
	}
#endif

#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
	DELETE_DESCRIPTOR(g_glslComputeController,
	                  OsdGLSLComputeController);
#endif

#undef DELETE_DESCRIPTOR
}
