#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vars.h"
#include "voronoi.h"

void initialize_qhull() {
  FILE *errorFile = stderr;

  QHULL_LIB_CHECK
  qh_zero(qh, errorFile);
}

void free_qhull() {
  FILE *errFile = stderr;
  int curlong, totlong;

  qh_freeqhull(qh, !qh_ALL);
  qh_memfreeshort(qh, &curlong, &totlong);
  if (curlong || totlong) {
    fprintf(errFile, "qhull internal warning: did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
  }

  free(qhReplicatedParticles);
}

void set_replicated_particles() {
  double enlargedBoxCoordMax;
  double enlargedBoxCoordMin;
  double enlargedBoxFrac;
  int count = 0;

#if DIMENSION == 2
  enlargedBoxFrac = 1.0;
#elif DIMENSION == 3
  if (numOfParticlesEachDim < 10) {
    enlargedBoxFrac = 1.0;
  } else if (numOfParticlesEachDim <= 64) {
    enlargedBoxFrac = ENLARGED_BOX_PARAM_A * pow(numOfParticlesEachDim, ENLARGED_BOX_PARAM_B);
  } else {
    enlargedBoxFrac = 0.2;
  }
#endif

  enlargedBoxCoordMax = boxSize * (1.0+enlargedBoxFrac);
  enlargedBoxCoordMin = -enlargedBoxFrac * boxSize;

  numOfReplicatedParticles = numOfParticles * pow(3, dimension);
  replicatedParticles = (struct particleData *)malloc(numOfReplicatedParticles * sizeof(struct particleData));
  if (replicatedParticles == NULL) {
    printf("Error: fail to allocate memory for struct particleData replicatedParticles.\n");
    exit(1);
  }

  for (int index = 0; index < numOfParticles; index++) {
    replicatedParticles[count++] = particles[index];
  }

#if DIMENSION == 2
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      if (i==0 && j==0) {
        continue;
      } else {
        for (int index = 0; index < numOfParticles; index++) {
          replicatedParticles[count] = particles[index];

          replicatedParticles[count].id = count;
          replicatedParticles[count].pos[0] += i*boxSize;
          replicatedParticles[count].pos[1] += j*boxSize;
          count++;
        }
      }
    }
  }
#elif DIMENSION == 3
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        if (i==0 && j==0 && k==0) {
          continue;
        } else {
          for (int index = 0; index < numOfParticles; index++) {
            replicatedParticles[count] = particles[index];

            replicatedParticles[count].id = count;
            replicatedParticles[count].pos[0] += i*boxSize;
            replicatedParticles[count].pos[1] += j*boxSize;
            replicatedParticles[count].pos[2] += k*boxSize;
            count++;
          }
        }
      }
    }
  }
#endif

  numOfSelectedReplicatedParticles = numOfReplicatedParticles;
  for (int i = 0; i < numOfReplicatedParticles; i++) {
    replicatedParticles[i].flagSelected = 1;
    for (int j = 0; j < dimension; j++) {
      if (replicatedParticles[i].pos[j]<enlargedBoxCoordMin || replicatedParticles[i].pos[j]>enlargedBoxCoordMax) {
        replicatedParticles[i].flagSelected = 0;
        numOfSelectedReplicatedParticles--;
        break;
      }
    }
  }
}

void replicated_particles_delaunay() {
  char flags[500];
  FILE *errorFile = stderr;
  FILE *outFile = NULL;
  int idArrayCount = 0, qhArrayCount = 0;
  int exitCode;
  boolT isMalloc = False;

  qhReplicatedParticles = (coordT *)malloc(numOfSelectedReplicatedParticles * dimension * sizeof(coordT));
  if (qhReplicatedParticles == NULL) {
    printf("Error: fail to allocate memory for array coordT qhReplicatedParticles.\n");
    exit(1);
  }
  selectedParticleId = (int *)malloc(numOfSelectedReplicatedParticles * sizeof(int));
  if (selectedParticleId == NULL) {
    printf("Error: fail to allocate memory for array int selectedParticleId.\n");
    exit(1);
  }

  for (int i = 0; i < numOfReplicatedParticles; i++) {
    if (replicatedParticles[i].flagSelected) {
      selectedParticleId[idArrayCount++] = i;
      for (int j = 0; j < dimension; j++) {
        qhReplicatedParticles[qhArrayCount++] = replicatedParticles[i].pos[j];
      }
    }
  }

  sprintf(flags, "qhull d QJ");
  exitCode = qh_new_qhull(qh, dimension, numOfSelectedReplicatedParticles, qhReplicatedParticles, 
                         isMalloc, flags, outFile, errorFile);
  if (exitCode) {
    printf("Error: fail to perform qhull.\n");
    exit(1);
  }

  qh_vertexneighbors(qh);
}

void map_replicated_particles_and_qhvertices() {
  int idTmp;
  int vertexIdMax;
  vertexT *vertex;

  vertexIdMax = 0;
  FORALLvertices {
    if (vertex->id > vertexIdMax) {
      vertexIdMax = vertex->id;
    }
  }

  qhVertexIdToReplicatedParticleId = (int *)malloc((vertexIdMax + 1) * sizeof(int));
  if (qhVertexIdToReplicatedParticleId == NULL) {
    printf("Error: fail to allocate memory for array qhVertexIdToReplicatedParticleId.\n");
    exit(1);
  }
  for (int i = 0; i < vertexIdMax; i++) {
    qhVertexIdToReplicatedParticleId[i] = -1;
  }

  replicatedParticleIdToQhVertex = malloc(numOfReplicatedParticles * sizeof(vertexT *));
  if (replicatedParticleIdToQhVertex == NULL) {
    printf("Error: fail to allocate memory for array replicatedParticleIdToQhVertex.\n");
    exit(1);
  }

  FORALLvertices {
    idTmp = qh_pointid(qh, vertex->point);
    idTmp = selectedParticleId[idTmp];
    qhVertexIdToReplicatedParticleId[vertex->id] = idTmp;
    replicatedParticleIdToQhVertex[idTmp] = vertex;
  }
}

void set_replicated_particle_neighbors() {
  int index;
  int neighborVertex_i, neighborVertex_n;
  unsigned int visitId;
  facetT *neighbor, **neighborp;
  vertexT *vertex;
  vertexT *neighborVertex;

  for (int i = 0; i < numOfReplicatedParticles; i++) {
    replicatedParticles[i].numOfNeighborParticles = 0;

    if (replicatedParticles[i].flagSelected) {
      visitId = ++qh->visit_id;
      vertex = replicatedParticleIdToQhVertex[i];
      vertex->visitid = visitId;
      FOREACHneighbor_(vertex) {
        if (neighbor->good) {
          FOREACHsetelement_i_(qh, vertexT, neighbor->vertices, neighborVertex) {
            if (neighborVertex->visitid != visitId) {
              neighborVertex->visitid = visitId;
              replicatedParticles[i].numOfNeighborParticles++;
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < numOfReplicatedParticles; i++) {
    if (replicatedParticles[i].flagSelected) {
      replicatedParticles[i].neighborParticleId = (int *)malloc(replicatedParticles[i].numOfNeighborParticles * \
                                                                sizeof(int));
      if (replicatedParticles[i].neighborParticleId == NULL) {
        printf("Error: fail to allocate memory for replicatedParticles[%d].neighborParticleId.\n", i);
        exit(1);
      }
      index = 0;

      visitId = ++qh->visit_id;
      vertex = replicatedParticleIdToQhVertex[i];
      vertex->visitid = visitId;
      FOREACHneighbor_(vertex) {
        if (neighbor->good) {
          FOREACHsetelement_i_(qh, vertexT, neighbor->vertices, neighborVertex) {
            if (neighborVertex->visitid != visitId) {
              neighborVertex->visitid = visitId;
              replicatedParticles[i].neighborParticleId[index++] = qhVertexIdToReplicatedParticleId[neighborVertex->id];
            }
          }
        }
      }
      sort_neighbor_particle_id(&replicatedParticles[i].neighborParticleId[0], replicatedParticles[i].numOfNeighborParticles);
    }
  }
}

void sort_neighbor_particle_id(int *ptr, int num) {
  int startIndex = 0;
  int endIndex = num - 1;
  int i = startIndex;
  int j = endIndex;
  int key = ptr[(i+j) / 2];
  int tmp;

  while (i < j) {
    for (; (i<endIndex) && (ptr[i]<key); i++);
    for (; (j>startIndex) && (ptr[j]>key); j--);
    if (i <= j) {
      tmp = ptr[i];
      ptr[i] = ptr[j];
      ptr[j] = tmp;
      i++;
      j--;
    }
  }

  if (i < endIndex) {
    sort_neighbor_particle_id(&ptr[i], endIndex - i + 1);
  }
  if (j > startIndex) {
    sort_neighbor_particle_id(&ptr[startIndex], j - startIndex + 1);
  }
}

void free_replicated_particle_neighbors() {
  for (int i = 0; i < numOfReplicatedParticles; i++) {
    if (replicatedParticles[i].flagSelected) {
      free(replicatedParticles[i].neighborParticleId);
    }
  }
}

void calculate_voronoi_volume_variance() {
  set_replicated_particles();
  replicated_particles_voronoi();
  calculate_each_voronoi_volume();
  calculate_volume_mean_and_variance();
  free_qhull();

  free(replicatedParticles);
  free(vertexPositions);
}

void replicated_particles_voronoi() {
  char flags[500];
  FILE *errFile = stderr;
  FILE *outFile = NULL;
  int qhArrayCount = 0;
  int exitCode;
  int numOfCenters;
  int numOfVertices = 1;
  boolT isLower;
  boolT isMalloc = False;
  facetT *facet;

  qhReplicatedParticles = (coordT *)malloc(numOfSelectedReplicatedParticles * dimension * sizeof(coordT));
  if (qhReplicatedParticles == NULL) {
    printf("Error: fail to allocate memory for array coordT qhReplicatedParticles.\n");
    exit(1);
  }
  for (int i = 0; i < numOfReplicatedParticles; i++) {
    if (replicatedParticles[i].flagSelected) {
      for (int j = 0; j < dimension; j++) {
        qhReplicatedParticles[qhArrayCount++] = replicatedParticles[i].pos[j];
      }
    }
  }

  qh_zero(qh, errFile); 
  sprintf(flags, "qhull s v QJ");  
  exitCode= qh_new_qhull(qh, dimension, numOfSelectedReplicatedParticles, qhReplicatedParticles, isMalloc, \
                      flags, outFile, errFile);
  if (exitCode) {
    printf("Error: fail to perform qhull in replicated_particles_voronoi().\n");
    exit(1);
  }

  qh_setvoronoi_all(qh);
  voronoiVertices = qh_markvoronoi(qh, qh->facet_list, NULL, !qh_ALL, &isLower, &numOfCenters);

  numOfFacets= (unsigned int) qh->num_facets;
  FORALLfacet_(qh->facet_list) {
    if (facet->visitid && facet->visitid < numOfFacets)
      numOfVertices++;
  }

  vertexPositions = (struct voronoiVertexPositionData *)malloc(numOfVertices * \
                    sizeof(struct voronoiVertexPositionData));
  if (vertexPositions == NULL) {
    printf("Error: fail to allocate momery for struct voronoiVertexPositionData vertexPositions.\n");
    exit(1);
  }

  for (int i = 0; i < dimension; i++)
    vertexPositions[0].pos[i] = qh_INFINITE;
 
  int i = 1;
  FORALLfacet_(qh->facet_list) {
    if (facet->visitid && facet->visitid < numOfFacets) {
      for (int j = 0; j < dimension; j++) {
        vertexPositions[i].pos[j] = facet->center[j];
      }
      i++;        
    }
  }
}

void calculate_each_voronoi_volume() {
  char flags[500];
  FILE *errFile = stderr;
  FILE *outFile = NULL;
  int curlong, totlong;
  int exitCode;
  int numOfNeighbors, numOfInfinityVertex;
  int particleId;
  int vertex_i, vertex_n;
  boolT isMalloc= False;
  coordT *voronoiCellVetex;
  facetT *facet;
  facetT *neighbor, **neighborp;
  qhT qh_qh_each_cell;
  qhT *qh_each_cell= &qh_qh_each_cell;
  vertexT *vertex;

  FOREACHvertex_i_(qh, voronoiVertices) {
    particleId = qh_pointid(qh, vertex->point);
    if (particleId < numOfParticles) {
      numOfNeighbors = 0;
      numOfInfinityVertex =0;

      if (vertex) {
        if (qh->hull_dim == 3)
          qh_order_vertexneighbors(qh, vertex);
        else if (qh->hull_dim >= 4)
          qsort(SETaddr_(vertex->neighbors, facetT), (size_t)qh_setsize(qh, vertex->neighbors), \
                sizeof(facetT *), qh_compare_facetvisit);
        FOREACHneighbor_(vertex) {
          if (neighbor->visitid == 0)
            numOfInfinityVertex = 1;
          else if (neighbor->visitid < numOfFacets)
            numOfNeighbors++;
        }
      }

      if (numOfInfinityVertex)
        numOfNeighbors++;

      voronoiCellVetex = (coordT *)malloc(numOfNeighbors * dimension * sizeof(coordT));
      if (voronoiCellVetex == NULL) {
        printf ("Error: fail to allocate memory for coordT voronoiCellVetex.\n");
        exit(1);
      }

      int i = 0;
      if (vertex) {
        FOREACHneighbor_(vertex) {
          if (neighbor->visitid == 0) {
            if (numOfInfinityVertex) {
              numOfInfinityVertex = 0;
              for (int j = 0; j < dimension; j++) {
                voronoiCellVetex[i] = vertexPositions[neighbor->visitid].pos[j];
                i++;
              } 
            }
          } else if (neighbor->visitid < numOfFacets) {
             for (int j = 0; j < dimension; j++) {
               voronoiCellVetex[i] = vertexPositions[neighbor->visitid].pos[j];
               i++;
             }
          }
        }
      }    

      qh_zero(qh_each_cell, errFile);  
      sprintf(flags, "qhull s");  
      exitCode = qh_new_qhull(qh_each_cell, dimension, numOfNeighbors, voronoiCellVetex, isMalloc,
                      flags, outFile, errFile);

      if (!exitCode) {
        facet = qh_each_cell->facet_list;
        qh_getarea(qh_each_cell, facet);
      }
      particles[particleId].volume = qh_each_cell->totvol;
      qh_freeqhull(qh_each_cell, !qh_ALL);
      qh_memfreeshort(qh_each_cell, &curlong, &totlong);
      if (curlong || totlong) {
        fprintf(errFile, "qhull internal warning: did not free %d bytes of long memory (%d pieces)\n", \
                totlong, curlong);
      }

      free(voronoiCellVetex);
    }
  }
}

void calculate_volume_mean_and_variance() {
  volumeMean = 0.0;
  volumeStd = 0.0;

  for (int i = 0; i < numOfParticles; i++) {
    volumeMean += particles[i].volume;
    volumeStd += particles[i].volume * particles[i].volume;
  }

  volumeMean /= numOfParticles;
  volumeStd = sqrt(volumeStd/numOfParticles - volumeMean*volumeMean);
}
