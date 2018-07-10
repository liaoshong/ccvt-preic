#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "libqhull_r/qhull_ra.h"

#include "ccdt.h"
#include "initialize.h"
#include "minimize.h"
#include "vars.h"
#include "voronoi.h"

void read_params(char *paramFilePath) {
  #define PARAM_MAX_NUM		10
  #define READ_MAX_NUM		300
  #define TYPE_INT		0
  #define TYPE_DOUBLE		1
  #define TYPE_STRING		2

  char str[READ_MAX_NUM], strTmp1[READ_MAX_NUM], strTmp2[READ_MAX_NUM];
  char paramName[PARAM_MAX_NUM][100];
  int paramCount = 0;
  int paramInputFlag[PARAM_MAX_NUM], paramType[PARAM_MAX_NUM];
  int pathFlag;
  void *paramAddr[PARAM_MAX_NUM];
  FILE *paramFile;

  strcpy(paramName[paramCount], "NUM_PART_EACH_DIM");
  paramType[paramCount] = TYPE_INT;
  paramAddr[paramCount] = &numOfParticlesEachDim;
  paramInputFlag[paramCount] = 0;
  paramCount++;

  strcpy(paramName[paramCount], "NUM_CAPACITY_EACH_DIM");
  paramType[paramCount] = TYPE_INT;
  paramAddr[paramCount] = &numOfCapacityEachDim;
  paramInputFlag[paramCount] = 0;
  paramCount++;

  strcpy(paramName[paramCount], "BOX_SIZE");
  paramType[paramCount] = TYPE_DOUBLE;
  paramAddr[paramCount] = &boxSize;
  paramInputFlag[paramCount] = 0;
  paramCount++;

  strcpy(paramName[paramCount], "RANDOM_SEED");
  paramType[paramCount] = TYPE_INT;
  paramAddr[paramCount] = &randomSeed;
  paramInputFlag[paramCount] = 0;
  paramCount++;

  strcpy(paramName[paramCount], "NUM_THREADS");
  paramType[paramCount] = TYPE_INT;
  paramAddr[paramCount] = &numOfThreads;
  paramInputFlag[paramCount] = 0;
  paramCount++;

  strcpy(paramName[paramCount], "OUTPUT_PATH");
  paramType[paramCount] = TYPE_STRING;
  paramAddr[paramCount] = outputPath;
  paramInputFlag[paramCount] = 0;
  paramCount++;

  paramFile = fopen(paramFilePath, "r");
  if (paramFile == NULL) {
    printf("Error: fail to open parameter file: %s.\n", paramFilePath);
    exit(1);
  }

  while (feof(paramFile) == 0) {
    fgets(str, READ_MAX_NUM, paramFile);
    if (str[0] == '%' || str[0] == 10) {
      continue;
    }
    int num = sscanf(str, "%s%s", strTmp1, strTmp2);
    if (num < 2) {
      continue;
    }

    for (int i = 0; i < paramCount; i++) {
      if (strcmp(strTmp1, paramName[i]) == 0) {
        if (paramInputFlag[i] == 0) {
          switch (paramType[i]) {
            case TYPE_INT:
              *((int *) paramAddr[i]) = atoi(strTmp2);
              break;
            case TYPE_DOUBLE:
              *((double *) paramAddr[i]) = atof(strTmp2);
              break;
            case TYPE_STRING:
              strcpy(paramAddr[i], strTmp2);
              break;
          }
          paramInputFlag[i] = 1;
        } else {
          printf("Error: %s in parameter file is re-defined.\n", strTmp1);
          exit(1);
        }
      }
    }
  }
  fclose(paramFile);

  for (int i = 0; i < paramCount; i++) {
    if (paramInputFlag[i] == 0) {
      printf("Error: %s is missing in parameter file.\n", paramName[i]);
      exit(1);
    }
  }

  dimension = DIMENSION;
  numOfParticles = pow(numOfParticlesEachDim, dimension);
  numOfCapacity = pow(numOfCapacityEachDim, dimension);
  numOfPoints = numOfCapacity * numOfParticles;
  numOfPointsEachDim = numOfCapacityEachDim * numOfParticlesEachDim;
  halfBoxSize = boxSize * 0.5;

  printf("INPUT PARAMETERS:\n");
  printf("  Number of particles on each dimension: %d\n", numOfParticlesEachDim);
  printf("  Capacity per particle per dimension: %d\n", numOfCapacityEachDim);
  printf("  Box size: %g\n", boxSize);
  printf("  Random seed: %d\n", randomSeed);
  printf("  Number of threads: %d\n", numOfThreads);
  printf("  Output path: %s\n\n", outputPath);

  if (DIMENSION <2 || DIMENSION >3) {
    printf("Dimension = %d\nError: currently, this code is for 2-D and 3-D cases only.\n", DIMENSION);
    exit(1);
  }
  if (numOfParticlesEachDim <= 1) {
    printf("Input NUM_PART_EACH_DIM = %d\nError: NUM_PART_EACH_DIM must be larger than 1.\n", numOfParticlesEachDim);
    exit(1);
  }
  if (numOfCapacityEachDim <= 1) {
    printf("Input NUM_CAPACITY_EACH_DIM = %d\nError: NUM_CAPACITY_EACH_DIM must be larger than 1.\n", \
           numOfCapacityEachDim);
    exit(1);
  }
  if (boxSize <= 0.0) {
    printf("Input BOX_SIZE = %g\nError: BOX_SIZE must be larger than 0.\n", boxSize);
    exit(1);
  }
  pathFlag = access(outputPath, F_OK);
  if (pathFlag) {
    printf("Error: given output path %s does not exist.\n", outputPath);
    exit(1);
  }
}

void initialize() {
  initialize_arrays();
#if CCDT_INIT == 1
  ccdt_initialize();
#endif
  assign_points_to_particles();
  initialize_bounding_radius();
}

void initialize_arrays() {
  double pointSeparation;
  int count = 0;

  pointSeparation = boxSize / numOfPointsEachDim;

  points = (struct pointData *)malloc(numOfPoints * sizeof(struct pointData));
  if (points == NULL) {
    printf("Error: fail to allocated memory for struct pointData points.\n");
    exit(1);
  }
#if DIMENSION == 2
  for (int i = 0; i < numOfPointsEachDim; i++) {
    int iMod2 = i % 2;
    for (int j = 0; j < numOfPointsEachDim; j++) {
      points[count].pos[0] = i * pointSeparation;
      points[count].pos[1] = (iMod2 ? numOfPointsEachDim-1-j : j) * pointSeparation;
      points[count].id = count;
      points[count].assignedParticleId = -1;
      count++;
    }
  }
#elif DIMENSION == 3
  for (int i = 0; i < numOfPointsEachDim; i++) {
    int iMod2 = i % 2;
    for (int j = 0; j < numOfPointsEachDim; j++) {
      int jMod2 = j % 2;
      for (int k = 0; k < numOfPointsEachDim; k++) {
        points[count].pos[0] = i * pointSeparation;
        points[count].pos[1] = (iMod2 ? numOfPointsEachDim-1-j : j) * pointSeparation;
        points[count].pos[2] = (jMod2 ? numOfPointsEachDim-1-k : k) * pointSeparation;
        points[count].id = count;
        points[count].assignedParticleId = -1;
        count++;
      }
    }
  }
#endif

  particles = (struct particleData *)malloc(numOfParticles * sizeof(struct particleData));
  if (particles == NULL) {
    printf("Error: fail to allocated memory for struct particleData particles.\n");
    exit(1);
  }
  srand(randomSeed);
  for (int i = 0; i < numOfParticles; i++) {
    particles[i].id = i;
    particles[i].numOfSwapsThisIter = 1;
    particles[i].capacity = 0;

    for (int j = 0; j < dimension; j++) {
      particles[i].pos[j] = (double) rand() / RAND_MAX * boxSize;
    }

    particles[i].pointsOfThisParticle = (struct pointList *)malloc(numOfCapacity * sizeof(struct pointList));
    if (particles[i].pointsOfThisParticle == NULL) {
      printf("Error: fail to allocate memory for struct pointList particles[%d].pointsOfThisParticle.\n", i);
      exit(1);
    }
  }

  numOfUnfilledParticles = numOfParticles;
  unfilledParticles = (int *)malloc(numOfUnfilledParticles * sizeof(int));

  for (int i = 0; i < numOfParticles; i++) {
    unfilledParticles[i] = i;
  }
}

void assign_points_to_particles() {
  int assignedParticleId;
  int idTmp;
  int particleId;
  boolT isOutside;
  coordT pointTmp[DIMENSION + 1];
  facetT *facet;
  realT bestDist, dist;
  vertexT *vertex, **vertexp;

  initialize_qhull();
  set_replicated_particles();
  replicated_particles_delaunay();
  map_replicated_particles_and_qhvertices();
  set_replicated_particle_neighbors();

  /* Assign the first point to its closest particle by using qh_findbestfacet() */
  for (int j = 0; j < dimension; j++) {
    pointTmp[j] = points[0].pos[j];
  }
  qh_setdelaunay(qh, dimension+1, 1, pointTmp);
  facet = qh_findbestfacet(qh, pointTmp, qh_ALL, &bestDist, &isOutside);
  bestDist = REALmax;
  FOREACHvertex_(facet->vertices) {
    dist = qh_pointdist(pointTmp, vertex->point, dimension);
    if (dist < bestDist) {
      bestDist = dist;
      points[0].initClosestParticleId = qhVertexIdToReplicatedParticleId[vertex->id];
    }
  }
  assignedParticleId = get_assigned_particle_id(points[0].initClosestParticleId);
  points[0].assignedParticleId = assignedParticleId;
  particles[assignedParticleId].pointsOfThisParticle[particles[assignedParticleId].capacity].pointId = 0;
  particles[assignedParticleId].capacity++;

  /* Assign other points to particles */
  for (int i = 1; i < numOfPoints; i++) {
    particleId = points[i-1].initClosestParticleId;
    bestDist = qh_pointdist(&points[i].pos[0], &replicatedParticles[particleId].pos[0], -dimension);
    points[i].initClosestParticleId = particleId;

    for (int j = 0; j < replicatedParticles[particleId].numOfNeighborParticles; j++) {
      idTmp = replicatedParticles[particleId].neighborParticleId[j];
      dist = qh_pointdist(&points[i].pos[0], &replicatedParticles[idTmp].pos[0], -dimension);
      if (dist < bestDist) {
        bestDist = dist;
        points[i].initClosestParticleId = idTmp;
      }
    }

    if (points[i].initClosestParticleId == points[i-1].initClosestParticleId) {
      assignedParticleId = get_assigned_particle_id(points[i-1].assignedParticleId);
    } else {
      assignedParticleId = get_assigned_particle_id(points[i].initClosestParticleId);
    }

    points[i].assignedParticleId = assignedParticleId;
    particles[assignedParticleId].pointsOfThisParticle[particles[assignedParticleId].capacity].pointId = i;
    particles[assignedParticleId].capacity++;

    if (particles[assignedParticleId].capacity == numOfCapacity) {
      int index = look_for_index(particles[assignedParticleId].id, 0, numOfUnfilledParticles-1);
      for (int j = index; j < numOfUnfilledParticles-1; j++) {
        unfilledParticles[j] = unfilledParticles[j+1];
      }
      numOfUnfilledParticles--;
    }
  }

  free(unfilledParticles);
  free(selectedParticleId);
  free(qhVertexIdToReplicatedParticleId);
  free(replicatedParticleIdToQhVertex);
  free_replicated_particle_neighbors();
  free(replicatedParticles);
  free_qhull();
}

int get_assigned_particle_id(int bestParticleId) {
  int idTmp = bestParticleId % numOfParticles;

  if (particles[idTmp].capacity < numOfCapacity) {
    return particles[idTmp].id;
  }

  for (int i = 0; i < replicatedParticles[bestParticleId].numOfNeighborParticles; i++) {
    idTmp = replicatedParticles[bestParticleId].neighborParticleId[i] % numOfParticles;

    if (particles[idTmp].capacity < numOfCapacity) {
      return particles[idTmp].id;
    }
  }

  return unfilledParticles[0];
}

int look_for_index(int id, int start, int end) {
  int middle;

  if(start > end)
    return -1;

  middle = (start + end) / 2;
  if(id == unfilledParticles[middle])
    return middle;
  else if(id < unfilledParticles[middle])
    return look_for_index(id, start, middle-1);
  else
    return look_for_index(id, middle+1, end);
}

void initialize_bounding_radius() {
  for (int i = 0; i < numOfParticles; i++) {
    calculate_particle_radius_and_energy(i);
  }
}
