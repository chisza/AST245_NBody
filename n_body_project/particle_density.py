import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math


def data_extraction(file_data):
	"""Extract the data into different arrays

	@param file_data: Text file loaded with numpy
	@type file_data: np.ndarray
	@return: tuple with the extracted data as arrays
	@rtype: tuple

	"""

	# particle number
	particle_number = file_data[:, 0]

	# particle total_mass
	mass = file_data[:, 1]

	# coordinates
	x_cord, y_cord, z_cord = file_data[:, 2], file_data[:, 3], file_data[:, 4]

	# velocities
	vx_vel, vy_vel, vz_vel = file_data[:, 5], file_data[:, 6], file_data[:, 7]

	# softening
	softening = file_data[:, 8]

	# potential
	potential = file_data[:, 9]

	print("Finished Extracting the Data")

	return particle_number, mass, x_cord, y_cord, z_cord, vx_vel, vy_vel, vz_vel, softening, potential


def hernquist_rho(ind_data, a):
	"""Calculate the density according to Hernquist

	@param ind_data: Tuple with array of radius and float of total total_mass
	@type ind_data: tuple[np.ndarray, np.ndarray]
	@param a: scale length
	@type a: float
	@return: density as an array
	@rtype: np.ndarray
	"""

	radius, tot_mass = ind_data

	rho = (tot_mass / (2 * np.pi)) * (a / radius) * (1 / ((radius + a) ** 3))

	return rho


def logathmic_shell_edges(x_cord, y_cord, z_cord, bin_number):
	"""Split the data into logarithmic shells / bins based on the radius"""

	radius = np.sqrt(x_cord ** 2 + y_cord **2 + z_cord ** 2)

	# calculate the start and end exponent, as the start value in logspace
	# is base ** value_given
	start_exp = math.floor(math.log10(np.min(radius)))
	end_exp = math.ceil(math.log10(np.max(radius)))

	# get logarithmically distributed bins
	bins = np.logspace(start_exp, end_exp, bin_number, base=10.0)
	widths = (bins[1:] - bins[:-1])

	return bins, radius, widths


def splitting_data_into_bins(bins, radius):
	"""Split the data into the shells"""

	particles_per_bin, bin_edges = np.histogram(radius, bins=bins)

	# check that the bin edges are still the same
	#for bin_edge in range(len(bin_edges)):
		#print(bin_edges[bin_edge], bins[bin_edge])

	return particles_per_bin, bin_edges


def masses_per_bin(particles_per_bin, masses):
	"""Calculate the mass enclosed in each bin"""

	# all particles have the same mass -> just multiply the number of
	# particles in each bin with the mass of a single particle
	masses_in_bin = particles_per_bin * np.min(masses)

	return masses_in_bin


def volume_per_bin(bin_edges):
	"""Calculate the volume for each bin as they are actually shells"""

	volume_in_bin = 4. / 3. * np.pi * (bin_edges[1:] ** 3 - bin_edges[:-1] ** 3)

	return volume_in_bin


def density_in_bin(masses_in_bin, volume_in_bin):
	"""Calculate the density in each bin"""

	density_in_bin = masses_in_bin / volume_in_bin

	return density_in_bin


def error_of_density(density_in_bin, shell_volumes, number_of_particles, masses):
	"""Calculate the error of the density distribution"""

	abs_error_of_density = (np.min(masses) / shell_volumes) * np.sqrt(number_of_particles)

	relative_error_of_density = np.divide(abs_error_of_density, density_in_bin, out=np.zeros_like(density_in_bin),
										  where=density_in_bin != 0)

	return abs_error_of_density, relative_error_of_density


def hernquist_fit(bin_edges, masses, density_in_bin, sigma = None):
	"""Fit the data to the Hernquist profile to get an estimate for a,
	then calculate the expected densities with appling the calculated a"""

	# get the average radius for each bin
	bin_radii = (bin_edges[1:] + bin_edges[:-1]) / 2

	tot_mass_system = np.empty(len(bin_radii))
	tot_mass_system.fill(np.sum(masses))
	#print(len(bin_radii))
	#print(len(density_in_bin))
	# get an estimate for scale length a
	inital_guess = 0.1

	param, param_cov = curve_fit(hernquist_rho, (bin_radii, tot_mass_system),
								 density_in_bin, p0=inital_guess, sigma=sigma, maxfev=100000)
	#print(param, param[0])

	# get the expected values based on the Hernquist profile and the estimated parameter
	expected_density_per_bin = hernquist_rho((bin_radii, tot_mass_system), param[0])

	return expected_density_per_bin, bin_radii


def error_measured_vs_expected(measured_density, expected_density_per_bin, particle_number, masses_per_shell, shell_volume):
	"""Compare the measured and the expected density of the particle distribution"""

	abs_error = np.mean(measured_density - expected_density_per_bin)

	rel_error = np.divide(abs_error, measured_density, out=np.zeros_like(measured_density), where= measured_density != 0)

	# expected number of particles
	exp_number_part = np.divide((expected_density_per_bin * shell_volume), masses_per_shell,
								out=np.zeros_like(masses_per_shell), where= masses_per_shell != 0)

	lambd = np.sqrt(np.mean(exp_number_part) / particle_number)

	return abs_error, rel_error, lambd


def plot_density_profile(x_cord, y_cord, z_cord, bin_number, masses):
	"""Plot the density measured based on the particles associated with the radii
	and compare it with the expected density based on the Hernquist profile"""

	# get the shell edges
	bins, radius, widths = logathmic_shell_edges(x_cord, y_cord, z_cord, bin_number)

	# split data into shells
	number_of_particles_per_bin, bin_edges = splitting_data_into_bins(bins, radius)

	# calculate the mass enclosed in each shell
	shell_masses = masses_per_bin(number_of_particles_per_bin, masses)

	# calculate the volume of each shell
	shell_volumes = volume_per_bin(bin_edges)

	# calculate the density in each shell
	shell_density_measured = density_in_bin(shell_masses, shell_volumes)

	# Calculate the error of the observed distribution
	abs_error, rel_error = error_of_density(shell_density_measured, shell_volumes, len(x_cord), masses)

	# calculate the density expected by the Hernquist profile
	expected_density_per_bin, bin_radii = hernquist_fit(bin_edges, masses, shell_density_measured, sigma=abs_error)

	# plot everything

	fig, ax = plt.subplots()

	ax.bar(bin_radii, shell_density_measured, width=widths,label="Number Density Measured")
	ax.errorbar(bin_radii, shell_density_measured, xerr=None, yerr=abs_error, color="orange", ls="none")
	ax.plot(bin_radii, expected_density_per_bin, c="red", label="Expected Number Density (Hernquist)")

	ax.set_title("Radial Number Density Profile")
	ax.legend()
	ax.set_xlabel("Radius")
	ax.set_ylabel("Number Density")
	ax.set_xscale("log")
	ax.set_yscale("log")

	print("Finished the density plot")

	return fig

