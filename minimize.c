#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "minimize.h"
#include "vars.h"

int minimize_energy() {
  struct groupData {
    int groupId;
    int numOfParticlesInThisGroup;
    int *particleIndex;
  };
  struct groupPairData {
    int group1;
    int group2;
  };

  struct groupData *groups;
  struct groupPairData *groupPairs;
  int addToGroupIndex;
  int count = 0;
  int numOfLoops = 2 * numOfThreads - 1;
  int numOfGroups = 2 * numOfThreads;
  int numOfGroupPairs = numOfThreads * (2*numOfThreads - 1);
  int numOfTotalSwaps = 0;
  int numOfAllocation = ceil((double)numOfParticles / (double)numOfGroups);

  groups = (struct groupData *)malloc(numOfGroups * sizeof(struct groupData));
  if (groups == NULL) {
    printf("Error: fail to allocate memory for struct groupData groups.\n");
    exit(1);
  }
  groupPairs = (struct groupPairData *)malloc(numOfGroupPairs * sizeof(struct groupPairData));
  if (groupPairs == NULL) {
    printf("Error: fail to allocate memory for struct groupPairData groupPairs.\n");
    exit(1);
  }

  sort_particle_number_of_swaps(&particles[0], numOfParticles);

  for (int i = 0; i < numOfGroups ; i++) {
    groups[i].groupId = i;
    groups[i].numOfParticlesInThisGroup = 0;
    groups[i].particleIndex = (int *)malloc(numOfAllocation * sizeof(int));
    if (groups[i].particleIndex == NULL) {
      printf("Error: fail to allocate memory for groups[%d].particleIndex.\n", i);
      exit(1);
    }
  }

  for (int i = 0; i < numOfParticles; i++) {
    addToGroupIndex = i % numOfGroups;
    groups[addToGroupIndex].particleIndex[groups[addToGroupIndex].numOfParticlesInThisGroup++] = i;
    particles[i].numOfSwapsLastIter = particles[i].numOfSwapsThisIter;
    particles[i].numOfSwapsThisIter = 0;
  }

  for (int i = 0; i < numOfLoops; i++) {
    for (int j = 0; j < numOfThreads; j++) {
      if (j == 0) {
        groupPairs[count].group1 = 2 * numOfThreads - 1;
        groupPairs[count].group2 = i;
      } else {
        groupPairs[count].group1 = (j + i) % (2*numOfThreads - 1);
        groupPairs[count].group2 = (2*numOfThreads - 1 - j + i) % (2*numOfThreads - 1);
      }
      count++;
    }
  }

  /* Particles in the same group */
  #pragma omp parallel for num_threads(numOfThreads)
  for (int groupIndex = 0; groupIndex < numOfGroups; groupIndex++) {
    double deltaEnergyTmp;
    double minEnergy1, minEnergy2;
    double particleSeparation;
    int index1, index2;
    int maxSwapNum;
    int numOfCandidates1, numOfCandidates2;
    int numOfSwapsThisGroup;
    int numOfSwapsThisThread = 0;
    int pointIdTmp, pointIdTmp1, pointIdTmp2;
    struct deltaEnergyData *deltaEnergyList1, *deltaEnergyList2;

    deltaEnergyList1 = (struct deltaEnergyData *)malloc(numOfCapacity * sizeof(struct deltaEnergyData));
    if (deltaEnergyList1 == NULL) {
      printf("Error: fail to allocate memory for struct deltaEnergyData deltaEnergyList1.\n");
      exit(1);
    }
    deltaEnergyList2 = (struct deltaEnergyData *)malloc(numOfCapacity * sizeof(struct deltaEnergyData));
    if (deltaEnergyList2 == NULL) {
      printf("Error: fail to allocate memory for struct deltaEnergyData deltaEnergyList2.\n");
      exit(1);
    }

    for (int i = 0; i < groups[groupIndex].numOfParticlesInThisGroup; i++) {
      index1 = groups[groupIndex].particleIndex[i];
      for (int j = i+1; j < groups[groupIndex].numOfParticlesInThisGroup; j++) {
        index2 = groups[groupIndex].particleIndex[j];
        particleSeparation = distance(&particles[index1].pos[0], &particles[index2].pos[0]);
        if ((particles[index1].numOfSwapsLastIter==0 && particles[index1].numOfSwapsThisIter==0 && \
            particles[index2].numOfSwapsLastIter==0 && particles[index2].numOfSwapsThisIter==0) || \
            particles[index1].boundingRadius + particles[index2].boundingRadius <= particleSeparation) {
          continue;
        }

        minEnergy1 = pow(particleSeparation, 2.0) - 2.0*particleSeparation*particles[index2].boundingRadius;
        numOfCandidates1 = 0;
        for (int k = 0; k < numOfCapacity; k++) {
          pointIdTmp = particles[index1].pointsOfThisParticle[k].pointId;
          deltaEnergyTmp = distance_square(&points[pointIdTmp].pos[0], &particles[index1].pos[0]) - \
                           distance_square(&points[pointIdTmp].pos[0], &particles[index2].pos[0]);
          if (deltaEnergyTmp > minEnergy1) {
            deltaEnergyList1[numOfCandidates1].pointIndex = k;
            deltaEnergyList1[numOfCandidates1].deltaEnergy = deltaEnergyTmp;
            numOfCandidates1++;
          }
        }

        if (numOfCandidates1 == 0) {
          continue;
        }
        sort_delta_energy(deltaEnergyList1, numOfCandidates1);

        minEnergy2 = -deltaEnergyList1[0].deltaEnergy;
        numOfCandidates2 = 0;
        for (int k = 0; k < numOfCapacity; k++) {
          pointIdTmp = particles[index2].pointsOfThisParticle[k].pointId;
          deltaEnergyTmp = distance_square(&points[pointIdTmp].pos[0], &particles[index2].pos[0]) - \
                           distance_square(&points[pointIdTmp].pos[0], &particles[index1].pos[0]);
          if (deltaEnergyTmp > minEnergy2) {
            deltaEnergyList2[numOfCandidates2].pointIndex = k;
            deltaEnergyList2[numOfCandidates2].deltaEnergy = deltaEnergyTmp;
            numOfCandidates2++;
          }
        }

        if (numOfCandidates2 == 0) {
          continue;
        }
        sort_delta_energy(deltaEnergyList2, numOfCandidates2);

        maxSwapNum = (numOfCandidates1 < numOfCandidates2) ? numOfCandidates1 : numOfCandidates2;
        numOfSwapsThisGroup = 0;
        for (int k = 0; k < maxSwapNum; k++) {
          if (deltaEnergyList1[k].deltaEnergy + deltaEnergyList2[k].deltaEnergy > 0.0) {
            pointIdTmp1 = particles[index1].pointsOfThisParticle[deltaEnergyList1[k].pointIndex].pointId;
            pointIdTmp2 = particles[index2].pointsOfThisParticle[deltaEnergyList2[k].pointIndex].pointId;

            particles[index1].pointsOfThisParticle[deltaEnergyList1[k].pointIndex].pointId = pointIdTmp2;
            particles[index1].pointsOfThisParticle[deltaEnergyList1[k].pointIndex].distanceSquare = \
              distance_square(&particles[index1].pos[0], &points[pointIdTmp2].pos[0]);
            particles[index2].pointsOfThisParticle[deltaEnergyList2[k].pointIndex].pointId = pointIdTmp1;
            particles[index2].pointsOfThisParticle[deltaEnergyList2[k].pointIndex].distanceSquare = \
              distance_square(&particles[index2].pos[0], &points[pointIdTmp1].pos[0]);

            particles[index1].numOfSwapsThisIter++;
            particles[index2].numOfSwapsThisIter++;
            numOfSwapsThisGroup++;
          } else {
            break;
          }
        }

        if (numOfSwapsThisGroup) {
          calculate_particle_radius_and_energy(index1);
          calculate_particle_radius_and_energy(index2);
          numOfSwapsThisThread += numOfSwapsThisGroup;
        }
      }
    }

    free(deltaEnergyList1);
    free(deltaEnergyList2);

    #pragma omp critical
    numOfTotalSwaps += numOfSwapsThisThread;
  }

  /* Particles in different groups */
  for (int loopIndex = 0; loopIndex < numOfLoops; loopIndex++) {
    #pragma omp parallel for num_threads(numOfThreads)
    for (int groupPairThisLoop = 0; groupPairThisLoop < numOfThreads; groupPairThisLoop++) {
      double deltaEnergyTmp;
      double minEnergy1, minEnergy2;
      double particleSeparation;
      int index1, index2;
      int maxSwapNum;
      int numOfCandidates1, numOfCandidates2;
      int numOfSwapsThisGroup;
      int numOfSwapsThisThread = 0;
      int pointIdTmp, pointIdTmp1, pointIdTmp2;
      int groupPairIndex = numOfThreads*loopIndex + groupPairThisLoop;
      int group1 = groupPairs[groupPairIndex].group1;
      int group2 = groupPairs[groupPairIndex].group2;
      struct deltaEnergyData *deltaEnergyList1, *deltaEnergyList2;

      deltaEnergyList1 = (struct deltaEnergyData *)malloc(numOfCapacity * sizeof(struct deltaEnergyData));
      if (deltaEnergyList1 == NULL) {
        printf("Error: fail to allocate memory for struct deltaEnergyData deltaEnergyList1.\n");
        exit(1);
      }
      deltaEnergyList2 = (struct deltaEnergyData *)malloc(numOfCapacity * sizeof(struct deltaEnergyData));
      if (deltaEnergyList2 == NULL) {
        printf("Error: fail to allocate memory for struct deltaEnergyData deltaEnergyList2.\n");
        exit(1);
      }

      for (int i = 0; i < groups[group1].numOfParticlesInThisGroup; i++) {
        index1 = groups[group1].particleIndex[i];
        for (int j = 0; j < groups[group2].numOfParticlesInThisGroup; j++) {
          index2 = groups[group2].particleIndex[j];
          particleSeparation = distance(&particles[index1].pos[0], &particles[index2].pos[0]);
          if ((particles[index1].numOfSwapsLastIter==0 && particles[index1].numOfSwapsThisIter==0 && \
              particles[index2].numOfSwapsLastIter==0 && particles[index2].numOfSwapsThisIter==0) || \
              particles[index1].boundingRadius + particles[index2].boundingRadius <= particleSeparation) {
            continue;
          }

          minEnergy1 = pow(particleSeparation, 2.0) - 2.0*particleSeparation*particles[index2].boundingRadius;
          numOfCandidates1 = 0;
          for (int k = 0; k < numOfCapacity; k++) {
            pointIdTmp = particles[index1].pointsOfThisParticle[k].pointId;
            deltaEnergyTmp = distance_square(&points[pointIdTmp].pos[0], &particles[index1].pos[0]) - \
                             distance_square(&points[pointIdTmp].pos[0], &particles[index2].pos[0]);
            if (deltaEnergyTmp > minEnergy1) {
              deltaEnergyList1[numOfCandidates1].pointIndex = k;
              deltaEnergyList1[numOfCandidates1].deltaEnergy = deltaEnergyTmp;
              numOfCandidates1++;
            }
          }

          if (numOfCandidates1 == 0) {
            continue;
          }
          sort_delta_energy(deltaEnergyList1, numOfCandidates1);

          minEnergy2 = -deltaEnergyList1[0].deltaEnergy;
          numOfCandidates2 = 0;
          for (int k = 0; k < numOfCapacity; k++) {
            pointIdTmp = particles[index2].pointsOfThisParticle[k].pointId;
            deltaEnergyTmp = distance_square(&points[pointIdTmp].pos[0], &particles[index2].pos[0]) - \
                             distance_square(&points[pointIdTmp].pos[0], &particles[index1].pos[0]);
            if (deltaEnergyTmp > minEnergy2) {
              deltaEnergyList2[numOfCandidates2].pointIndex = k;
              deltaEnergyList2[numOfCandidates2].deltaEnergy = deltaEnergyTmp;
              numOfCandidates2++;
            }
          }

          if (numOfCandidates2 == 0) {
            continue;
          }
          sort_delta_energy(deltaEnergyList2, numOfCandidates2);

          maxSwapNum = (numOfCandidates1 < numOfCandidates2) ? numOfCandidates1 : numOfCandidates2;
          numOfSwapsThisGroup = 0;
          for (int k = 0; k < maxSwapNum; k++) {
            if (deltaEnergyList1[k].deltaEnergy + deltaEnergyList2[k].deltaEnergy > 0.0) {
              pointIdTmp1 = particles[index1].pointsOfThisParticle[deltaEnergyList1[k].pointIndex].pointId;
              pointIdTmp2 = particles[index2].pointsOfThisParticle[deltaEnergyList2[k].pointIndex].pointId;

              particles[index1].pointsOfThisParticle[deltaEnergyList1[k].pointIndex].pointId = pointIdTmp2;
              particles[index1].pointsOfThisParticle[deltaEnergyList1[k].pointIndex].distanceSquare = \
                distance_square(&particles[index1].pos[0], &points[pointIdTmp2].pos[0]);
              particles[index2].pointsOfThisParticle[deltaEnergyList2[k].pointIndex].pointId = pointIdTmp1;
              particles[index2].pointsOfThisParticle[deltaEnergyList2[k].pointIndex].distanceSquare = \
                distance_square(&particles[index2].pos[0], &points[pointIdTmp1].pos[0]);

              particles[index1].numOfSwapsThisIter++;
              particles[index2].numOfSwapsThisIter++;
              numOfSwapsThisGroup++;
            } else {
              break;
            }
          }

          if (numOfSwapsThisGroup) {
            calculate_particle_radius_and_energy(index1);
            calculate_particle_radius_and_energy(index2);
            numOfSwapsThisThread += numOfSwapsThisGroup;
          }
        }
      }

      free(deltaEnergyList1);
      free(deltaEnergyList2);

      #pragma omp critical
      numOfTotalSwaps += numOfSwapsThisThread;
    }
  }

  #pragma omp parallel for num_threads(numOfThreads)
  for(int i = 0; i < numOfParticles; i++) {
    move_particle_to_centroid(i);
    calculate_particle_radius_and_energy(i);
  }

  free(groups);
  free(groupPairs);

  return numOfTotalSwaps;
}

void sort_particle_number_of_swaps(struct particleData *ptr, int num) {
  int startIndex = 0;
  int endIndex = num - 1;
  int i = startIndex;
  int j = endIndex;
  int key = ptr[(i+j) / 2].numOfSwapsThisIter;
  struct particleData temp;

  while (i < j) {
    for (; (i<endIndex) && (ptr[i].numOfSwapsThisIter<key); i++);
    for (; (j>startIndex) && (ptr[j].numOfSwapsThisIter>key); j--);
    if (i <= j) {
      temp = ptr[i];
      ptr[i] = ptr[j];
      ptr[j] = temp;
      i++;
      j--;
    }
  }

  if (i < endIndex) {
    sort_particle_number_of_swaps(&ptr[i], endIndex - i + 1);
  }
  if (j > startIndex) {
    sort_particle_number_of_swaps(&ptr[startIndex], j - startIndex + 1);
  }
}

void sort_delta_energy(struct deltaEnergyData *ptr, int num) {
  int startIndex = 0;
  int endIndex = num - 1;
  int i = startIndex;
  int j = endIndex;
  double key = ptr[(i+j) / 2].deltaEnergy;
  struct deltaEnergyData temp;

  while (i < j) {
    for (; (i<endIndex) && (ptr[i].deltaEnergy>key); i++);
    for (; (j>startIndex) && (ptr[j].deltaEnergy<key); j--);
    if (i <= j) {
      temp = ptr[i];
      ptr[i] = ptr[j];
      ptr[j] = temp;
      i++;
      j--;
    }
  }

  if (i < endIndex) {
    sort_delta_energy(&ptr[i], endIndex - i + 1);
  }
  if (j > startIndex) {
    sort_delta_energy(&ptr[startIndex], j - startIndex + 1);
  }
}

void move_particle_to_centroid(int particleId) {
  double firstPointPos[DIMENSION];
  int firstPointId = particles[particleId].pointsOfThisParticle[0].pointId;

  for (int i = 0; i < dimension; i++) {
    particles[particleId].pos[i] = 0.0;
    firstPointPos[i] = points[firstPointId].pos[i];
  }
  for (int i = 0; i < numOfCapacity; i++) {
    for (int j = 0; j < dimension; j++) {
      int pointId = particles[particleId].pointsOfThisParticle[i].pointId;
      particles[particleId].pos[j] += periodic_pos(firstPointPos[j], points[pointId].pos[j]);
    }
  }
  for (int i = 0; i < dimension; i++) {
    particles[particleId].pos[i] = periodic_pos(halfBoxSize, particles[particleId].pos[i]/numOfCapacity);
  }
}

double periodic_pos(double x0, double x) {
  double dx = x - x0;

  if (dx < -halfBoxSize) {
    return x + boxSize;
  } else if (dx > halfBoxSize) {
    return x - boxSize;
  } else {
    return x;
  }
}

void calculate_particle_radius_and_energy(int particleId) {
  double distanceSquare;

  particles[particleId].boundingRadius = -1.0;
  particles[particleId].energy = 0.0;
  for (int i = 0; i < particles[particleId].capacity; i++) {
    distanceSquare = distance_square(&particles[particleId].pos[0], \
                       &points[particles[particleId].pointsOfThisParticle[i].pointId].pos[0]);
    particles[particleId].pointsOfThisParticle[i].distanceSquare = distanceSquare;
    particles[particleId].energy += distanceSquare;
    if (distanceSquare > particles[particleId].boundingRadius) {
      particles[particleId].boundingRadius = distanceSquare;
    }
  }

  particles[particleId].boundingRadius = sqrt(particles[particleId].boundingRadius);
}

double distance(double *pos1, double *pos2) {
  return sqrt(distance_square(pos1, pos2));
}

double distance_square(double *pos1, double *pos2) {
  double distanceSquare = 0.0;
  double dx;

  for (int i = 0; i < dimension; i++) {
    dx = periodic_dx(pos1[i], pos2[i]);
    distanceSquare += dx * dx;
  }

  return distanceSquare;
}

double periodic_dx(double x1, double x2) {
  double dx = x1 - x2;

  if (dx < -halfBoxSize) {
    return dx + boxSize;
  } else if (dx > halfBoxSize) {
    return dx - boxSize;
  } else {
    return dx;
  }
}

double calculate_total_energy() {
  double totalEnergy = 0.0;

  for (int i = 0; i < numOfParticles; i++) {
    totalEnergy += particles[i].energy;
  }

  return totalEnergy;
}
