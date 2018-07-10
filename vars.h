#ifndef VARS_H
#define VARS_H

#include "libqhull_r/qhull_ra.h"
#include "configure.h"

struct deltaEnergyData {
  int pointIndex;
  double deltaEnergy;
};

struct neighborFacetData {
  int id;
  int vertexAId;
  int vertexBId;
  int vertexCId;
  double vertexAPos[DIMENSION];
  double vertexBPos[DIMENSION];
  double vertexCPos[DIMENSION];
  double paramA;
  double paramB;
  double paramC;
  double paramD;
};

struct pointList {
  int pointId;
  double distanceSquare;
};

struct particleData {
  int id;
  int numOfSwapsThisIter;
  int numOfSwapsLastIter;
  double pos[DIMENSION];
  double boundingRadius;
  double energy;
  double voronoiVolume;
  int capacity;
  struct pointList *pointsOfThisParticle;
  int numOfNeighborParticles;
  int *neighborParticleId;
  int numOfNeighborFacets;
  struct neighborFacetData *neighborFacets;
  unsigned visitId;
  double volume;
  int flagSelected;
};

struct pointData {
  int id;
  double pos[DIMENSION];
  int assignedParticleId;
  int initClosestParticleId;
};

struct voronoiVertexPositionData {
  double pos[DIMENSION];
};

extern struct particleData *particles, *replicatedParticles;
extern struct pointData *points;
extern struct voronoiVertexPositionData *vertexPositions;

extern char outputPath[1000];

extern double boxSize, halfBoxSize;
extern double volumeMean, volumeStd;

extern int dimension;
extern int numOfCapacityEachDim, numOfCapacity;
extern int numOfFacets;
extern int numOfParticlesEachDim, numOfParticles;
extern int numOfPointsEachDim, numOfPoints;
extern int numOfReplicatedParticles;
extern int numOfSelectedReplicatedParticles;
extern int numOfThreads;
extern int numOfUnfilledParticles;
extern int randomSeed;
extern int *qhVertexIdToReplicatedParticleId;
extern int *selectedParticleId;
extern int *unfilledParticles;

extern coordT *qhReplicatedParticles;
extern qhT qh_qh;
extern qhT *qh;
extern setT *voronoiVertices;
extern vertexT **replicatedParticleIdToQhVertex;

#endif
