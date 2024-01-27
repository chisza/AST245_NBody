import math
import matplotlib.axes
import matplotlib.figure

from numba import jit, prange
import numpy as np
import matplotlib.pyplot as plt


@jit(nopython=True, parallel=True)
def direct_force_calculation(mass, x_cord, y_cord, z_cord, softening):
	"""Calculate the forces with a brute force calculation

	@param mass:  total_mass of the particles as an array
	@type mass: np.ndarray
	@param x_cord: x coordinates as an array
	@type x_cord: np.ndarray
	@param y_cord: y coordinates as an array
	@type y_cord: np.ndarray
	@param z_cord: z coordinates as an array
	@type z_cord: np.ndarray
	@param softening: softening parameter as array
	@type softening: np.ndarray
	@return: tuple of absolute forces, absolute radii, and tuple of accelerations
	@rtype: tuple[np.ndarray, np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]

	"""

	G = 1.0
	N = len(x_cord)
	# initialize empty arrays for the acceleration
	ax_accel, ay_accel, az_accel = np.zeros_like(x_cord), np.zeros_like(y_cord), np.zeros_like(z_cord)

	# initialize empty arrays for the forces
	Fx, Fy, Fz, abs_F = np.zeros_like(x_cord), np.zeros_like(y_cord), np.zeros_like(z_cord), np.zeros_like(x_cord)

	if isinstance(softening, int):
		softening = [softening * len(x_cord)]

	print("Begin brute force calculation")
	print(softening)

	for i in prange(N):
		# per = i / N * 100
		# print(per)
		for j in range(N):
			if i == j:
				continue
			rx = x_cord[j] - x_cord[i]
			ry = y_cord[j] - y_cord[i]
			rz = z_cord[j] - z_cord[i]

			abs_r = np.sqrt(rx ** 2 + ry ** 2 + rz ** 2)

			ax_accel[i] += -G * (mass[j] / (abs_r ** 2 + softening[i] ** 2) ** (3 / 2)) * rx
			ay_accel[i] += -G * (mass[j] / (abs_r ** 2 + softening[i] ** 2) ** (3 / 2)) * ry
			az_accel[i] += -G * (mass[j] / (abs_r ** 2 + softening[i] ** 2) ** (3 / 2)) * rz

	Fx = ax_accel * mass
	Fy = ay_accel * mass
	Fz = az_accel * mass
	abs_F = np.sqrt(Fx**2 + Fy**2 + Fz**2)

	abs_r = np.sqrt(x_cord ** 2 + y_cord ** 2 + z_cord ** 2)

	accelerations = (ax_accel, ay_accel, az_accel)

	print(len(abs_F))

	print("Finished brute force calculation")

	return abs_F, abs_r, accelerations


@jit(nopython=True)
def half_mass_radius(x_coordinate, y_coordinate, z_coordinate, mass):
	"""Calculate the half-total_mass radius

	@param x_coordinate: x coordinates as an array
	@type x_coordinate: np.ndarray
	@param y_coordinate: y coordinates as an array
	@type y_coordinate: np.ndarray
	@param z_coordinate: z coordinates as an array
	@type z_coordinate: np.ndarray
	@param mass: total_mass of particles as an array
	@type mass: np.ndarray
	@return: Half total_mass radius
	@rtype: float

	"""

	print("Calculate the half total_mass radius")

	# get the radius for each particle and sort them
	radii = np.sqrt(x_coordinate**2 + y_coordinate**2 + z_coordinate**2)
	radii = np.sort(radii)

	# calculate the total total_mass of the system
	tot_mass = np.sum(mass)

	# calculate half of the total_mass
	half_mass = tot_mass / 2.
	acc_mass = 0
	# all particles have the same mass
	for i in range(len(x_coordinate)):
		if acc_mass <= half_mass:
			acc_mass += mass[i]
			ind = i
		else:
			break

	# get the radius corresponding with half the total_mass
	half_mass_r = radii[ind]
	half_mass = acc_mass

	return half_mass_r, half_mass


def particles_in_hmr(half_mass_radius, x_coordinate, y_coordinate, z_coordinate):
	"""Calculate the particles that are within the half total_mass radius

	@param half_mass_radius: The half total_mass radius of the particle distribution
	@type half_mass_radius: float
	@param x_coordinate: x coordinates of the particles as array
	@type x_coordinate: np.ndarray
	@param y_coordinate: y coordinates of the particles as array
	@type y_coordinate: np.ndarray
	@param z_coordinate: z coordinates of the particles as array
	@type z_coordinate: np.ndarray
	@return: Array of the positions of each particle (x, y, z)
	@rtype: np.ndarray

	"""

	print("Get the particles within the half total_mass radius")

	# get a list of the particles within the half total_mass radius
	# get that for all coordinates
	# is the complete radius, that is relevant, not the single coordinates
	# get the radius for each particle and sort them
	x_coordinate = np.sort(x_coordinate)
	y_coordinate = np.sort(y_coordinate)
	z_coordinate = np.sort(z_coordinate)

	radii = np.sqrt(x_coordinate ** 2 + y_coordinate ** 2 + z_coordinate ** 2)
	radii = np.sort(radii)

	# Replace the loop for filtering particles with array slicing
	mask = radii <= half_mass_radius
	positions = np.array([x_coordinate[mask], y_coordinate[mask], z_coordinate[mask]]).T

	return positions


@jit(nopython=True, parallel=True)
def mean_inter_particle_separation(positions):
	"""Calculate the mean inter particle separation

	@param positions: Positions of the particles as array (x, y, z)
	@type positions: np.ndarray
	@return: the mean inter particle separation
	@rtype: float

	"""

	print("Calculate the mean inter particle distance")

	# calculate from the remaining particles the mean inter particle distance
	number_part = len(positions)
	sum_dist_arr = np.zeros_like(positions)

	for i in prange(number_part):
		#per = i / number_part * 100
		#print(per)
		for j in range(number_part):
			if i == j:
				continue
			rx = positions[j][0] - positions[i][0]
			ry = positions[j][1] - positions[i][1]
			rz = positions[j][2] - positions[i][2]

			distance = np.sqrt(rx**2 + ry**2 + rz**2)
			#print(distance)
			sum_dist_arr[i] += distance

	sum_dist = np.sum(sum_dist_arr)
	mean_dist = sum_dist / (number_part * (number_part - 1))

	return mean_dist


def softening_values(mean_dist, exp, step_size):
	"""Get a list of different softening values, based on the mean inter particle separation

	@param mean_dist: the mean inter particle separation
	@type mean_dist: float
	@param exp: the exponent, resulting in a list
	@type exp: int
	@param step_size: in the exponents, which difference there is
	@type step_size: int
	@return: list of softening values
	@rtype: list
	"""

	exp_of_mean_dist = math.floor(math.log10(mean_dist))
	n = 10 ** exp_of_mean_dist
	exp = list(range(-exp, exp + 1, step_size))
	softs = []
	for i in exp:
		softs.append(n * (10 ** i))

	return softs


def softening_plot(soft_values, particle_number, mean_int, mass, x_cord, y_cord, z_cord, softening):
	"""Get a plot for different softenings

	@param soft_values: list of softenings values
	@type soft_values: list
	@param particle_number: Particle indexes
	@type particle_number: np.ndarray
	@param mean_int: mean inter particle separation
	@type mean_int: float
	@param mass: particle total_mass as array
	@type mass: np.ndarray
	@param x_cord: x coordinates of particles as array
	@type x_cord: np.ndarray
	@param y_cord: y coordinates of particles as array
	@type y_cord: np.ndarray
	@param z_cord: z coordinates of particles as array
	@type z_cord: np.ndarray
	@param softening: original softening values as array
	@type softening: np.ndarray
	@return: tuple with the figure of the softening plot, and a tuple with the plot and axis for the radius plot
	@rtype: tuple[matplotlib.figure.Figure, tuple[matplotlib.figure.Figure, plt.Axes]]
	"""

	print("Start softening plot")
	# Start plot
	soft_plot, ax_soft = plt.subplots()
	rad_plot, ax_rad = plt.subplots()

	# get a list of the lenght of the particles, so that only every 100th particle is shown
	par = list(range(0, len(particle_number), 1))

	for i in soft_values:
		softening[:] = i
		F_abs, r_abs, accel = direct_force_calculation(mass, x_cord, y_cord, z_cord, softening)
		F_abs_sort = np.rec.fromarrays([F_abs, r_abs], dtype=np.dtype([('F_abs', np.float32), ('r_abs', np.float32)])) # find the permutation to sort r_abs fro lowest to highest
		F_abs_sort.sort()

		r_abs_sort = F_abs_sort.r_abs
		F_abs_sort = F_abs_sort.F_abs

		if i / (10 ** math.floor(math.log10(mean_int))) == 1:
			ax_soft.scatter(par[::100], F_abs_sort[::100], color="red", edgecolors='black', label=f'Softening={i:.0e}')
			ax_rad.scatter(r_abs_sort[::100], F_abs_sort[::100], color="red", edgecolors='black', label=f'Softening={i:.0e}')
		else:
			ax_soft.scatter(par[::100], F_abs_sort[::100], label=f'Softening={i:.0e}')
			ax_rad.scatter(r_abs_sort[::100], F_abs_sort[::100], label=f'Softening={i:.0e}')

	ax_soft.legend()
	#plt.xscale("log")
	ax_soft.set_yscale("log")
	ax_soft.set_xlabel('Particle Index')
	ax_soft.set_ylabel('Order of Force')

	ax_rad.legend()
	ax_rad.set_xscale("log")
	ax_rad.set_yscale("log")
	ax_rad.set_xlabel('absolute radius')
	ax_rad.set_ylabel('Order of Force')

	print("Finished softening plot")

	return soft_plot, (rad_plot, ax_rad)


