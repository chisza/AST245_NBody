# do not use numba here, it does NOT work
import numpy as np
import math
import sys
import matplotlib.pyplot as plt

def center_of_mass(x_coordinates, y_coordinates, z_coordinates, masses):
	"""Calculate the center of total_mass coordinates of a particle distribution

	@param x_coordinates: x coordinates as array
	@type x_coordinates: np.ndarray
	@param y_coordinates: y coordinates as array
	@type y_coordinates: np.ndarray
	@param z_coordinates: z coordinates as array
	@type z_coordinates: np.ndarray
	@param masses: masses as array
	@type masses: np.ndarray
	@return: center of mass as tuple
	@rtype: tuple
	"""

	# calculate the total total_mass of the system
	tot_mass = np.sum(masses)

	r_x = np.sum(masses * x_coordinates) / tot_mass
	r_y = np.sum(masses * y_coordinates) / tot_mass
	r_z = np.sum(masses * z_coordinates) / tot_mass

	return (r_x, r_y, r_z)


def max_radius(x_coordinates, y_coordinates, z_coordinates):
	"""Find the overall radius of the system

	@param x_coordinates: x coordinates as array
	@type x_coordinates: np.ndarray
	@param y_coordinates: y coordinates as array
	@type y_coordinates: np.ndarray
	@param z_coordinates: z coordinates as array
	@type z_coordinates: np.ndarray
	@return: tuple of the absolute radii of the particles and the overall absolute radius of the system
	@rtype: tuple
	"""

	radii = np.sqrt(x_coordinates ** 2 + y_coordinates ** 2 + z_coordinates ** 2)

	maximal_radius = max(radii)

	return radii, maximal_radius


def split_data_into_shells(x_coordinates, y_coordinates, z_coordinates, n_shells):
	"""Split the data into different shells

	@param x_coordinates: x coordinates as array
	@type x_coordinates: np.ndarray
	@param y_coordinates: y coordinates as array
	@type y_coordinates: np.ndarray
	@param z_coordinates: z coordinates as array
	@type z_coordinates: np.ndarray
	@param n_shells: the number of shells to be generated
	@type n_shells: int
	@return: tuple of the index of particles in shell, the shell boundaries, the shell thickness
	@rtype: tuple
	"""

	# Calculate the maximal radius
	radii, maximal_radius = max_radius(x_coordinates, y_coordinates, z_coordinates)

	start_exp = math.floor(math.log10(np.min(radii)))
	end_exp = math.ceil(math.log10(np.max(radii)))

	# get logarithmically distributed shells
	# add + 1 to generate the number of shells specified
	# this needs n + 1 boundaries
	shell_boundaries = np.logspace(start_exp, end_exp, n_shells + 1, base=10.0)

	# get an array with the shell thicknesses of every shell
	# the shell thicknesses also have a lograrithmic spacing
	# with this procedure, every two boundaries are subtracted from each other
	# get N (50) thicknesses for 51 boundaries generated for 50 shells
	shell_thickness = shell_boundaries[1:] - shell_boundaries[:-1]


	in_shells = []

	for i in range(n_shells):
		# Find indices of particles within the current shell
		lower_bound = shell_boundaries[i]
		upper_bound = shell_boundaries[i + 1]
		in_shell = np.where((radii >= lower_bound) & (radii <= upper_bound))
		in_shells.append(in_shell[0])  # Use [0] to get the actual indices

	# get the maximal length of the arrays with the shell radii
	max_length = max([len(ele) for ele in in_shells])

	# prepare an array that contains -1 with the max_length to fill it
	# so that all arrays have the same length in the end
	in_shells_new = np.empty((n_shells, max_length))
	in_shells_new.fill(-1) # the index of a particle is never negative

	# fill the new array
	for i, ele in enumerate(in_shells):
		in_shells_new[i, :len(ele)] = ele

	return in_shells_new, shell_boundaries, shell_thickness


def find_shell_index(r, shell_boundaries):
	# This function finds the shell index for a given radius 'r'
	# It associates a radius with a shell
	# It returns the shell a particle is in
	for i in range(len(shell_boundaries) - 1):
		if shell_boundaries[i] <= r < shell_boundaries[i + 1]:
			return i
	return -1  # Return -1 if the radius is outside all shells


def rho_of_shell(shell_indices, shell_boundaries, part_ind, masses):
	"""Find the density of a shell, calculate the total_mass of each shell
	and the volume of each shell, and the density is then rho = total_mass / V
	per shell
	"""

	# shell_indices is an array, that contains arrays(= shells)
	# with the indices of the particles in that shell

	# Find the indices of particles in each shell using boolean indexing

	mass_of_shells = []
	count = 0

	# check if a particle is in a shell
	# if yes, then add the mass of the particle to the total mass of that shell
	for shell in range(np.shape(shell_indices)[0]):
		mass_shell = 0
		shell = np.array(shell_indices[shell], dtype=np.int64)
		for particle in shell:
			if particle > -1:
				count += 1
				mass_shell += masses[particle]
		mass_of_shells.append(mass_shell)

	#print(mass_of_shells)
	#print(count)

	shell_volumes = []

	# calculate the shell volumes
	# these should be logarithmic, as the boundaries have been calculated logarithmically
	# QUESTION as the boundaries are logarithmic, the volumes should be logarithmic as well, right?
	for i in range(len(shell_boundaries) - 1):
		lower_boundary = shell_boundaries[i]
		upper_boundary = shell_boundaries[i + 1]
		shell_volume = 4. / 3. * np.pi * (upper_boundary ** 3 - lower_boundary ** 3)
		shell_volumes.append(shell_volume)

	#print(shell_volumes)

	mass_of_shells = np.array(mass_of_shells)
	shell_volumes = np.array(shell_volumes)
	rho_of_shells = mass_of_shells / shell_volumes
	#print(len(rho_of_shells))

	return rho_of_shells



def shell_masses(rho_of_shells, shell_boundaries, shell_thickness):
	"""Calculate the total_mass of each shell"""

	masses_of_spheres = []

	shell_boundaries = shell_boundaries[:]

	# calculate the density for each shell
	# QUESTION even when the innermost shell boundary is not 0 anymore because of logspace, can I still sum up from the first and not 0th element?
	# I would say yes, the innermost boundary is the starting point, needed to create shells but not a radius
	# I need to sum up towards to
	shell_density_masses = shell_thickness * shell_boundaries[1:] ** 2 * rho_of_shells

	# QUESTION this includes the mass density of the shell the particle is in
	# is that how it works?
	# If the shell the particle is in is not included, what happens
	# with the particles in the inner most shell?
	# Or have just many more shells to be chosen?

	# the number of boundaries for the shells is one bigger than
	# the number of shells -> to loop over just the shells -> -1
	for i in range(len(shell_boundaries)-1):
		# sum the densities up to the desired shell according to the
		# formula in galactic dynamics
		summe = 4 * np.pi * np.sum(shell_density_masses[:i+1])
		masses_of_spheres.append(summe)

	# calculate the masses enclosed in each sphere bounded by a shell_boundary
	#for i in range(len(shell_boundaries)-1): # for 50 shells if 51 boundaries, comes from another function
		#the_summed_part = np.sum(shell_thickness * shell_boundaries[1:i+1] ** 2 * rho_of_shells[:i])
		#total_mass_of_system = 4 * np.pi * the_summed_part
		#masses_of_spheres.append(total_mass_of_system)

	#print(masses_of_spheres)
	#print(len(masses_of_spheres))

	masses_of_spheres = np.array(masses_of_spheres)

	return masses_of_spheres


def forces_when_looking_at_shells(shell_masses, shell_boundaries, x_coordinates, y_coordinates, z_coordinates):
	"""Calculate the forces on a particle when subdividing the system into shells
	and looking at the inner shell for the force"""

	G = 1

	# gives the absolute radii for each particle
	radii, maximal_radius = max_radius(x_coordinates, y_coordinates, z_coordinates)

	# calculate the unit vectors for all particles

	# first combine the x, y, and z coordinates for each particle
	particle_positions = np.stack((x_coordinates, y_coordinates, z_coordinates), axis=-1)

	# calculate the unit vector for each particle
	unit_vector = particle_positions / radii[:, np.newaxis]

	absolute_forces = []

	# get the according inner shell masses for each particle
	# check of the particle radius is smaller than the shell radius
	# the first shell boundary is 0 because it is the center,
	# so it is not really a boundary, so it has to be excluded
	# to retrieve the effective radius for each shell
	shell_boundaries = shell_boundaries[1:]

	count = 0
	for particle in range(len(x_coordinates)):
		#print(particle)
		# this shell index, it corresponds to the shell the particle is in
		# to calculate the force a particle experiences -> shell_ind - 1
		# to get the inner shell
		shell_ind = np.max(np.where(radii[particle] <= shell_boundaries))
		#print(shell_ind)
		#if shell_ind >= 0:
			#count += 1
	#print(count)

		if shell_ind > 0:
			force = - (G * shell_masses[shell_ind]) / radii[particle] ** 2 * unit_vector[particle]
			abs_force = np.linalg.norm(force)
			#print(abs_force)
			absolute_forces.append(abs_force)
		else:
			absolute_forces.append(0)

	absolute_forces = np.array(absolute_forces)
	# print("here")
	# print(len(absolute_forces))
	# print(absolute_forces)
	# print(radii)

	return absolute_forces, radii










