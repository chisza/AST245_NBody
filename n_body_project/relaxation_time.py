import numpy as np
from numba import jit


@jit(nopython=True)
def typical_particle_velocity(half_mass, half_mass_radius):
	"""calculate the typcial particle velocity based on the
	half mass radius"""

	print("Calculating typical particle velocity")

	G = 1

	typical_velocity = np.sqrt((G * half_mass) / half_mass_radius)

	return typical_velocity


@jit(nopython=True)
def crossing_time(half_mass_radius, typical_velocity):
	"""Calculate the crossing time of the particle
	based on the half mass radius and the typical particle velocity"""

	print("Calculating crossing time")

	crossing_time = half_mass_radius / typical_velocity

	return crossing_time


@jit(nopython=True)
def relaxation_time(number_of_particles_half_mass, half_mass, half_mass_radius):
	"""Calculate the relaxation time of the system based on the half mass radius"""

	print("Calculating relaxation time")

	typical_velocity = typical_particle_velocity(half_mass, half_mass_radius)

	crossing_time_system = crossing_time(half_mass_radius, typical_velocity)

	# the numpy function log is the natural logarithm
	relaxation_time_system = number_of_particles_half_mass / (8 * np.log(number_of_particles_half_mass)) * crossing_time_system

	return relaxation_time_system, crossing_time_system
