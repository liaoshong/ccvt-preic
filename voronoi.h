#ifndef VORONOI_H
#define VORONOI_H

#define ENLARGED_BOX_PARAM_A		4.2
#define ENLARGED_BOX_PARAM_B		-0.772

void calculate_each_voronoi_volume();
void calculate_volume_mean_and_variance();
void calculate_voronoi_volume_variance();
void free_qhull();
void free_replicated_particle_neighbors();
void initialize_qhull();
void map_replicated_particles_and_qhvertices();
void replicated_particles_delaunay();
void replicated_particles_voronoi();
void set_replicated_particles();
void set_replicated_particle_neighbors();
void sort_neighbor_particle_id(int *ptr, int num);

#endif
