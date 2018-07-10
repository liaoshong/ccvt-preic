#include <omp.h>
#include <stdio.h>
#include <string.h>
#include "initialize.h"
#include "minimize.h"
#include "output.h"
#include "vars.h"
#include "voronoi.h"

int main(int argc, char *argv[]) {
  char paramFilePath[1000];
  double minimizeEnergyWtimeStart, minimizeEnergyWtimeEnd, wtimeStart, wtimeEnd;
  int iterCount = 0;
  int numOfSwaps;

  if (argc < 2) {
    printf("Error: arguments missing.\n");
    printf("Specify the parameter file in the command line, e.g., ./ccvt-preic <parameterFile>\n");
    exit(1);
  } else {
    strcpy(paramFilePath, argv[1]);
  }

  wtimeStart = omp_get_wtime();
  read_params(paramFilePath);
  initialize();
  wtimeEnd = omp_get_wtime();
  printf("Initialize. Wall time used: %g seconds.\n\n", wtimeEnd - wtimeStart);

  minimizeEnergyWtimeStart = omp_get_wtime();
  do {
    printf("Iteration %d. ", ++iterCount);
    numOfSwaps = minimize_energy();
    printf("Swap number: %d, total energy: %g.\n", numOfSwaps, calculate_total_energy());
  } while(numOfSwaps);
  minimizeEnergyWtimeEnd = omp_get_wtime();
  printf("Optimize. Wall time used: %g seconds.\n\n", minimizeEnergyWtimeEnd - minimizeEnergyWtimeStart);

  calculate_voronoi_volume_variance();
  printf("Compute Voronoi volume. Mean: %g, Std: %g (%.2f%%).\n\n", volumeMean, volumeStd, \
         volumeStd/volumeMean*100.);

  output();

  wtimeEnd = omp_get_wtime();
  printf("Program ends. Total wall time used: %g seconds.\n", wtimeEnd - wtimeStart);

  return 0;
}
