#include <stdio.h>
#include <stdlib.h>

#include "output.h"
#include "vars.h"

void output() {
  char fileName[500];
  FILE *out;

  sprintf(fileName, "%s/ccvt_particle_%d_capacity_%d.txt", outputPath, numOfParticlesEachDim, numOfCapacityEachDim);
  out = fopen(fileName, "w");
  if (out == NULL) {
    printf("Error: fail to open file: %s.\n", fileName);
    exit(1);
  }

#if DIMENSION == 2
  fprintf(out, "x y\n");
#elif DIMENSION == 3
  fprintf(out, "x y z\n");
#endif
  for (int i = 0; i < numOfParticles; i++) {
    for (int j = 0; j < dimension; j++) {
      fprintf(out, "%g ", particles[i].pos[j]);
    }
    fprintf(out, "\n");
  }
  fclose(out);

#if DIMENSION == 3
  output_gadget_format();
#endif
}

void output_gadget_format() {
  struct gadgetHeader {
    int npart[6];
    double mass[6];
    double time;
    double redshift;
    int flagSfr;
    int flagFeedback;
    unsigned int npartTotal[6];
    int flagCooling;
    int numFiles;
    double boxSize;
    double omega0;
    double omegaLambda;
    double hubbleParam;
    int flagStellarage;
    int flagMetals;
    unsigned int npartTotalHighWord[6];
    int  flagEntropyInsteadU;
    char fill[60];
  };

  struct gadgetHeader header;
  char fileName[500];
  FILE *out;
  float tmp;
  int blockSize;
  unsigned idTmp;

  for (int i = 0; i < 6; i++) {
    header.npart[i] = 0;
    header.npartTotal[i] = 0;
    header.mass[i] = 0.0;
    header.npartTotalHighWord[i] = 0;
  }
  header.npart[1] = header.npartTotal[1] = numOfParticles;
  header.mass[1] = 27.7536 * pow(boxSize, 3.) / numOfParticles;

  header.time = 1.0;
  header.redshift = 0.0;
  header.flagSfr = 0;
  header.flagFeedback = 0;
  header.flagCooling = 0;
  header.numFiles = 1;
  header.boxSize = boxSize;
  header.omega0 = 1.0;
  header.omegaLambda = 0.0;
  header.hubbleParam = 0.70;
  header.flagStellarage = 0;
  header.flagMetals = 0;
  header.flagEntropyInsteadU = 0;

  sprintf(fileName, "%s/gadget_particle_%d_capacity_%d", outputPath, numOfParticlesEachDim, numOfCapacityEachDim);
  out = fopen(fileName, "wb");
  if (out == NULL) {
    printf("Error: fail to open file: %s.\n", fileName);
    exit(1);
  }

  blockSize = sizeof(header);
  fwrite(&blockSize, 1, sizeof(blockSize), out);
  fwrite(&header, 1, sizeof(header), out);
  fwrite(&blockSize, 1, sizeof(blockSize), out);

  blockSize = sizeof(float) * numOfParticles * 3;
  fwrite(&blockSize, 1, sizeof(blockSize), out);
  for (int i = 0; i < numOfParticles; i++) {
    for (int j = 0; j < 3; j++) {
      tmp = (float)particles[i].pos[j];
      fwrite(&tmp, 1, sizeof(float), out);
    }
  }
  fwrite(&blockSize, 1, sizeof(blockSize), out);

  tmp = 0.0;
  fwrite(&blockSize, 1, sizeof(blockSize), out);
  for (int i = 0; i < numOfParticles; i++) {
    for (int j = 0; j < 3; j++) {
      fwrite(&tmp, 1, sizeof(float), out);
    }
  }
  fwrite(&blockSize, 1, sizeof(blockSize), out);

  blockSize = sizeof(unsigned) * numOfParticles;
  fwrite(&blockSize, 1, sizeof(blockSize), out);
  for (int i = 0; i < numOfParticles; i++) {
    idTmp = i;
    fwrite(&idTmp, 1, sizeof(unsigned), out);
  }
  fwrite(&blockSize, 1, sizeof(blockSize), out);

  fclose(out);
}
