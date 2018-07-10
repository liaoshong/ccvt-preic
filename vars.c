#include "vars.h"

struct particleData *particles, *replicatedParticles;
struct pointData *points;
struct voronoiVertexPositionData *vertexPositions;

char outputPath[1000];

double boxSize, halfBoxSize;
double volumeMean, volumeStd;

int dimension;
int numOfCapacityEachDim, numOfCapacity;
int numOfFacets;
int numOfParticlesEachDim, numOfParticles;
int numOfPointsEachDim, numOfPoints;
int numOfReplicatedParticles;
int numOfSelectedReplicatedParticles;
int numOfThreads;
int numOfUnfilledParticles;
int randomSeed;
int *qhVertexIdToReplicatedParticleId;
int *selectedParticleId;
int *unfilledParticles;

coordT *qhReplicatedParticles;
qhT qh_qh;
qhT *qh = &qh_qh;
setT *voronoiVertices;
vertexT **replicatedParticleIdToQhVertex;
