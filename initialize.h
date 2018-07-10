#ifndef INITIALIZE_H
#define INITIALIZE_H

void assign_points_to_particles();
int get_assigned_particle_id(int bestParticleId);
void initialize_arrays();
void initialize();
void initialize_bounding_radius();
int look_for_index(int id, int start, int end);
void read_params(char *paramFilePath);

#endif
