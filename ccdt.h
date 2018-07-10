#ifndef CCDT_H
#define CCDT_H

#define NUM_DELAUNAY_ITER		2
#define NUM_VARIANCE_RELAX_ITER		30

void ccdt_initialize();
void free_replicated_particle_neighbor_facets();
void relax_area_or_volume_variance();
void set_replicated_particle_neighbor_facets();
void update_particle_pos();

#endif
