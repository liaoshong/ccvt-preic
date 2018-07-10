#ifndef MINIMIZE_H
#define MINIMIZE_H

#include "vars.h"

void calculate_particle_radius_and_energy(int particleId);
double calculate_total_energy();
double distance(double *pos1, double *pos2);
double distance_square(double *pos1, double *pos2);
int minimize_energy();
void move_particle_to_centroid(int particleId);
double periodic_dx(double x1, double x2);
double periodic_pos(double x0, double x);
void sort_delta_energy(struct deltaEnergyData *ptr, int num);
void sort_particle_number_of_swaps(struct particleData *ptr, int num);

#endif
