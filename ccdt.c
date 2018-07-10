#include <stdlib.h>
#include "libqhull_r/qhull_ra.h"

#include "ccdt.h"
#include "minimize.h"
#include "voronoi.h"

void ccdt_initialize() {
  int numOfDelaunayIter = NUM_DELAUNAY_ITER;

  for (int iter = 0; iter < numOfDelaunayIter; iter++) {
    initialize_qhull();
    set_replicated_particles();
    replicated_particles_delaunay();
    map_replicated_particles_and_qhvertices();
    set_replicated_particle_neighbor_facets();
    relax_area_or_volume_variance();
    update_particle_pos();

    free(selectedParticleId);
    free(qhVertexIdToReplicatedParticleId);
    free(replicatedParticleIdToQhVertex);
    free_replicated_particle_neighbor_facets();
    free(replicatedParticles);
    free_qhull();
  }
}

void set_replicated_particle_neighbor_facets() {
  double det;
  double vectors[DIMENSION][DIMENSION];
  double vertexPos[DIMENSION][DIMENSION];
  int facetIndex, vertexIdIndex;
  int neighborVertex_i, neighborVertex_n;
  int vertexId[DIMENSION];
  facetT *neighbor, **neighborp;
  vertexT *vertex;
  vertexT *neighborVertex;

  for (int i = 0; i < numOfReplicatedParticles; i++) {
    replicatedParticles[i].numOfNeighborFacets = 0;

    if (replicatedParticles[i].flagSelected) {
      vertex = replicatedParticleIdToQhVertex[i];
      FOREACHneighbor_(vertex) {
        if (neighbor->good) {
          replicatedParticles[i].numOfNeighborFacets++;
        }
      }

      replicatedParticles[i].neighborFacets = (struct neighborFacetData *)malloc(replicatedParticles[i].numOfNeighborFacets * \
                                               sizeof(struct neighborFacetData));
      if (replicatedParticles[i].neighborFacets == NULL) {
        printf("Error: fail to allocate memory for replicatedParticles[%d].neighborFacets.\n", i);
        exit(1);
      }

      facetIndex = 0;
      FOREACHneighbor_(vertex) {
        if (neighbor->good) {
          vertexIdIndex = 0;
          FOREACHsetelement_i_(qh, vertexT, neighbor->vertices, neighborVertex) {
            if (qhVertexIdToReplicatedParticleId[neighborVertex->id] != i) {
              vertexId[vertexIdIndex++] = qhVertexIdToReplicatedParticleId[neighborVertex->id];
            }
          }

          #if DIMENSION == 2
          vertexPos[0][0] = replicatedParticles[vertexId[0]].pos[0];
          vertexPos[0][1] = replicatedParticles[vertexId[0]].pos[1];
          vertexPos[1][0] = replicatedParticles[vertexId[1]].pos[0];
          vertexPos[1][1] = replicatedParticles[vertexId[1]].pos[1];

          vectors[0][0] = vertexPos[0][0] - replicatedParticles[i].pos[0];
          vectors[0][1] = vertexPos[0][1] - replicatedParticles[i].pos[1];
          vectors[1][0] = vertexPos[1][0] - replicatedParticles[i].pos[0];
          vectors[1][1] = vertexPos[1][1] - replicatedParticles[i].pos[1];

          det = vectors[0][0]*vectors[1][1] - vectors[0][1]*vectors[1][0];
          if (det > 0.0) {
            replicatedParticles[i].neighborFacets[facetIndex].vertexAId = vertexId[0];
            replicatedParticles[i].neighborFacets[facetIndex].vertexBId = vertexId[1];
          } else {
            replicatedParticles[i].neighborFacets[facetIndex].vertexAId = vertexId[1];
            replicatedParticles[i].neighborFacets[facetIndex].vertexBId = vertexId[0];
          }
          #elif DIMENSION == 3
          vertexPos[0][0] = replicatedParticles[vertexId[0]].pos[0];
          vertexPos[0][1] = replicatedParticles[vertexId[0]].pos[1];
          vertexPos[0][2] = replicatedParticles[vertexId[0]].pos[2];

          vertexPos[1][0] = replicatedParticles[vertexId[1]].pos[0];
          vertexPos[1][1] = replicatedParticles[vertexId[1]].pos[1];
          vertexPos[1][2] = replicatedParticles[vertexId[1]].pos[2];

          vertexPos[2][0] = replicatedParticles[vertexId[2]].pos[0];
          vertexPos[2][1] = replicatedParticles[vertexId[2]].pos[1];
          vertexPos[2][2] = replicatedParticles[vertexId[2]].pos[2];

          vectors[0][0] = replicatedParticles[i].pos[0] - vertexPos[0][0];
          vectors[0][1] = replicatedParticles[i].pos[1] - vertexPos[0][1];
          vectors[0][2] = replicatedParticles[i].pos[2] - vertexPos[0][2];

          vectors[1][0] = vertexPos[1][0] - vertexPos[0][0];
          vectors[1][1] = vertexPos[1][1] - vertexPos[0][1];
          vectors[1][2] = vertexPos[1][2] - vertexPos[0][2];

          vectors[2][0] = vertexPos[2][0] - vertexPos[0][0];
          vectors[2][1] = vertexPos[2][1] - vertexPos[0][1];
          vectors[2][2] = vertexPos[2][2] - vertexPos[0][2];

          det = vectors[0][0] * (vectors[1][1]*vectors[2][2] - vectors[2][1]*vectors[1][2]) + \
                vectors[0][1] * (vectors[1][2]*vectors[2][0] - vectors[2][2]*vectors[1][0]) + \
                vectors[0][2] * (vectors[1][0]*vectors[2][1] - vectors[2][0]*vectors[1][1]); 
          if (det > 0.0) {
            replicatedParticles[i].neighborFacets[facetIndex].vertexAId = vertexId[0];
            replicatedParticles[i].neighborFacets[facetIndex].vertexBId = vertexId[1];
            replicatedParticles[i].neighborFacets[facetIndex].vertexCId = vertexId[2];
          } else {
            replicatedParticles[i].neighborFacets[facetIndex].vertexAId = vertexId[0];
            replicatedParticles[i].neighborFacets[facetIndex].vertexBId = vertexId[2];
            replicatedParticles[i].neighborFacets[facetIndex].vertexCId = vertexId[1];
          }
          #endif

          facetIndex++;
        }
      }
    }
  }
}

void free_replicated_particle_neighbor_facets() {
  for (int i = 0; i < numOfReplicatedParticles; i++) {
    if (replicatedParticles[i].flagSelected) {
      free(replicatedParticles[i].neighborFacets);
    }
  }
}

void relax_area_or_volume_variance() {
  int numOfRelaxations = NUM_VARIANCE_RELAX_ITER;

  for (int iter = 0; iter < numOfRelaxations; iter++) {
    #if DIMENSION == 2
    double deltaX, deltaY, xNew, yNew;
    double numReciprocal;
    double sumAiAi, sumAiBi, sumAiCi, sumBiBi, sumBiCi, sumAmAn, sumAmBn, sumAmCn, sumBmBn, sumBmCn;
    double u11, u12, u21, u22, v1, v2;
    double xA, xB, yA, yB;
    int idA, idB;

    for (int i = 0; i < numOfParticles; i++) {
      for (int j = 0; j < replicatedParticles[i].numOfNeighborFacets; j++) {
        idA = replicatedParticles[i].neighborFacets[j].vertexAId;
        idB = replicatedParticles[i].neighborFacets[j].vertexBId;
        for (int k = 0; k < dimension; k++) {
          replicatedParticles[i].neighborFacets[j].vertexAPos[k] = replicatedParticles[idA].pos[k];
          replicatedParticles[i].neighborFacets[j].vertexBPos[k] = replicatedParticles[idB].pos[k];
        }
       
        xA = replicatedParticles[i].neighborFacets[j].vertexAPos[0];
        yA = replicatedParticles[i].neighborFacets[j].vertexAPos[1];
        xB = replicatedParticles[i].neighborFacets[j].vertexBPos[0];
        yB = replicatedParticles[i].neighborFacets[j].vertexBPos[1];

        replicatedParticles[i].neighborFacets[j].paramA = yA - yB;
        replicatedParticles[i].neighborFacets[j].paramB = xB - xA;
        replicatedParticles[i].neighborFacets[j].paramC = xA*yB - xB*yA;
      }

      sumAiAi = 0.0;
      sumAiBi = 0.0;
      sumAiCi = 0.0;
      sumBiBi = 0.0;
      sumBiCi = 0.0;
      sumAmAn = 0.0;
      sumAmBn = 0.0;
      sumAmCn = 0.0;
      sumBmBn = 0.0;
      sumBmCn = 0.0;

      for (int j = 0; j < replicatedParticles[i].numOfNeighborFacets; j++) {
        sumAiAi += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[j].paramA;
        sumAiBi += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[j].paramB;
        sumAiCi += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[j].paramC;
        sumBiBi += replicatedParticles[i].neighborFacets[j].paramB * replicatedParticles[i].neighborFacets[j].paramB;
        sumBiCi += replicatedParticles[i].neighborFacets[j].paramB * replicatedParticles[i].neighborFacets[j].paramC;

        for (int k = 0; k < replicatedParticles[i].numOfNeighborFacets; k++) {
          sumAmAn += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[k].paramA;
          sumAmBn += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[k].paramB;
          sumAmCn += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[k].paramC;
          sumBmBn += replicatedParticles[i].neighborFacets[j].paramB * replicatedParticles[i].neighborFacets[k].paramB;
          sumBmCn += replicatedParticles[i].neighborFacets[j].paramB * replicatedParticles[i].neighborFacets[k].paramC;
        }
      }

      numReciprocal = 1.0 / replicatedParticles[i].numOfNeighborFacets;

      u11 = sumAiAi - sumAmAn*numReciprocal;
      u12 = sumAiBi - sumAmBn*numReciprocal;
      u21 = u12;
      u22 = sumBiBi - sumBmBn*numReciprocal;
      v1 = sumAmCn*numReciprocal - sumAiCi;
      v2 = sumBmCn*numReciprocal - sumBiCi;

      xNew = (u22*v1 - u12*v2) / (u11*u22 - u12*u21);
      yNew = (u11*v2 - u21*v1) / (u11*u22 - u12*u21);

      deltaX = xNew - replicatedParticles[i].pos[0];
      deltaY = yNew - replicatedParticles[i].pos[1];

      for (int j = i; j < numOfReplicatedParticles; j += numOfParticles) {
        replicatedParticles[j].pos[0] += deltaX;
        replicatedParticles[j].pos[1] += deltaY;
      }
    }

    #elif DIMENSION == 3
    double A11, A12, A21, A22, B1, B2;
    double deltaX, deltaY, deltaZ, xNew, yNew, zNew;
    double numReciprocal;
    double sumAiAi, sumAiBi, sumAiCi, sumAiDi, sumBiBi, sumBiCi, sumBiDi, sumCiCi, sumCiDi, \
           sumAmAn, sumAmBn, sumAmCn, sumAmDn, sumBmBn, sumBmCn, sumBmDn, sumCmCn, sumCmDn;
    double u11, u12, u13, u21, u22, u23, u31, u32, u33, v1, v2, v3;
    double xA, xB, xC, yA, yB, yC, zA, zB, zC;
    int idA, idB, idC;

    for (int i = 0; i < numOfParticles; i++) {
      for (int j = 0; j < replicatedParticles[i].numOfNeighborFacets; j++) {
        idA = replicatedParticles[i].neighborFacets[j].vertexAId;
        idB = replicatedParticles[i].neighborFacets[j].vertexBId;
        idC = replicatedParticles[i].neighborFacets[j].vertexCId;
        for (int k = 0; k < dimension; k++) {
          replicatedParticles[i].neighborFacets[j].vertexAPos[k] = replicatedParticles[idA].pos[k];
          replicatedParticles[i].neighborFacets[j].vertexBPos[k] = replicatedParticles[idB].pos[k];
          replicatedParticles[i].neighborFacets[j].vertexCPos[k] = replicatedParticles[idC].pos[k];
        }
       
        xA = replicatedParticles[i].neighborFacets[j].vertexAPos[0];
        yA = replicatedParticles[i].neighborFacets[j].vertexAPos[1];
        zA = replicatedParticles[i].neighborFacets[j].vertexAPos[2];

        xB = replicatedParticles[i].neighborFacets[j].vertexBPos[0];
        yB = replicatedParticles[i].neighborFacets[j].vertexBPos[1];
        zB = replicatedParticles[i].neighborFacets[j].vertexBPos[2];

        xC = replicatedParticles[i].neighborFacets[j].vertexCPos[0];
        yC = replicatedParticles[i].neighborFacets[j].vertexCPos[1];
        zC = replicatedParticles[i].neighborFacets[j].vertexCPos[2];

        replicatedParticles[i].neighborFacets[j].paramA = (yB - yA)*(zC - zA) - (yC - yA)*(zB - zA);
        replicatedParticles[i].neighborFacets[j].paramB = (zB - zA)*(xC - xA) - (zC - zA)*(xB - xA);
        replicatedParticles[i].neighborFacets[j].paramC = (xB - xA)*(yC - yA) - (xC - xA)*(yB - yA);
        replicatedParticles[i].neighborFacets[j].paramD = -1.0 * (xA*replicatedParticles[i].neighborFacets[j].paramA + \
                                                                  yA*replicatedParticles[i].neighborFacets[j].paramB + \
                                                                  zA*replicatedParticles[i].neighborFacets[j].paramC);
      }

      sumAiAi = 0.0;
      sumAiBi = 0.0;
      sumAiCi = 0.0;
      sumAiDi = 0.0;
      sumBiBi = 0.0;
      sumBiCi = 0.0;
      sumBiDi = 0.0;
      sumCiCi = 0.0;
      sumCiDi = 0.0;

      sumAmAn = 0.0;
      sumAmBn = 0.0;
      sumAmCn = 0.0;
      sumAmDn = 0.0; 
      sumBmBn = 0.0;
      sumBmCn = 0.0;
      sumBmDn = 0.0;
      sumCmCn = 0.0;
      sumCmDn = 0.0;

      for (int j = 0; j < replicatedParticles[i].numOfNeighborFacets; j++) {
        sumAiAi += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[j].paramA;
        sumAiBi += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[j].paramB;
        sumAiCi += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[j].paramC;
        sumAiDi += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[j].paramD;
        sumBiBi += replicatedParticles[i].neighborFacets[j].paramB * replicatedParticles[i].neighborFacets[j].paramB;
        sumBiCi += replicatedParticles[i].neighborFacets[j].paramB * replicatedParticles[i].neighborFacets[j].paramC;
        sumBiDi += replicatedParticles[i].neighborFacets[j].paramB * replicatedParticles[i].neighborFacets[j].paramD;
        sumCiCi += replicatedParticles[i].neighborFacets[j].paramC * replicatedParticles[i].neighborFacets[j].paramC;
        sumCiDi += replicatedParticles[i].neighborFacets[j].paramC * replicatedParticles[i].neighborFacets[j].paramD;

        for (int k = 0; k < replicatedParticles[i].numOfNeighborFacets; k++) {
          sumAmAn += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[k].paramA;
          sumAmBn += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[k].paramB;
          sumAmCn += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[k].paramC;
          sumAmDn += replicatedParticles[i].neighborFacets[j].paramA * replicatedParticles[i].neighborFacets[k].paramD;
          sumBmBn += replicatedParticles[i].neighborFacets[j].paramB * replicatedParticles[i].neighborFacets[k].paramB;
          sumBmCn += replicatedParticles[i].neighborFacets[j].paramB * replicatedParticles[i].neighborFacets[k].paramC;
          sumBmDn += replicatedParticles[i].neighborFacets[j].paramB * replicatedParticles[i].neighborFacets[k].paramD;
          sumCmCn += replicatedParticles[i].neighborFacets[j].paramC * replicatedParticles[i].neighborFacets[k].paramC;
          sumCmDn += replicatedParticles[i].neighborFacets[j].paramC * replicatedParticles[i].neighborFacets[k].paramD;
        }
      }

      numReciprocal = 1.0 / replicatedParticles[i].numOfNeighborFacets;

      u11 = sumAiAi - sumAmAn*numReciprocal;
      u12 = sumAiBi - sumAmBn*numReciprocal;
      u13 = sumAiCi - sumAmCn*numReciprocal;
      u21 = u12;
      u22 = sumBiBi - sumBmBn*numReciprocal;
      u23 = sumBiCi - sumBmCn*numReciprocal;
      u31 = u13;
      u32 = u23;
      u33 = sumCiCi - sumCmCn*numReciprocal;

      v1 = sumAmDn*numReciprocal - sumAiDi;
      v2 = sumBmDn*numReciprocal - sumBiDi;
      v3 = sumCmDn*numReciprocal - sumCiDi;

      A11 = u11*u33 - u13*u31;
      A12 = u12*u33 - u13*u32;
      A21 = u21*u33 - u23*u31;
      A22 = u22*u33 - u23*u32;
      B1 = v1*u33 - v3*u13;
      B2 = v2*u33 - v3*u23;

      xNew = (A22*B1 - A12*B2) / (A11*A22 - A12*A21);
      yNew = (A11*B2 - A21*B1) / (A11*A22 - A12*A21);
      zNew = (v1 - u11*xNew - u12*yNew) / u13;

      deltaX = xNew - replicatedParticles[i].pos[0];
      deltaY = yNew - replicatedParticles[i].pos[1];
      deltaZ = zNew - replicatedParticles[i].pos[2];

      for (int j = i; j < numOfReplicatedParticles; j += numOfParticles) {
        replicatedParticles[j].pos[0] += deltaX;
        replicatedParticles[j].pos[1] += deltaY;
        replicatedParticles[j].pos[2] += deltaZ;
      }
    }
    #endif
  }
}

void update_particle_pos() {
  for (int i = 0; i < numOfParticles; i++) {
    for (int j = 0; j < dimension; j++) {
      if (replicatedParticles[i].pos[j] < 0.0) {
        particles[i].pos[j] = replicatedParticles[i].pos[j] + boxSize;
      } else if (replicatedParticles[i].pos[j] >= boxSize) {
        particles[i].pos[j] = replicatedParticles[i].pos[j] - boxSize;
      } else {
        particles[i].pos[j] = replicatedParticles[i].pos[j];
      }
    }
  }
}
