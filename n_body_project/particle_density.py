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


def particle_density(mass, x_coordinate, y_coordinate, z_coordinate, bin_number):
	"""Get a plot of the density of the particles and add the estimated
	density of Hernquist

	@param mass: Mass as an array
	@type mass: np.ndarray
	@param x_coordinate: x coordinates as an array
	@type x_coordinate: np.ndarray
	@param y_coordinate: y coordinates as an array
	@type y_coordinate: np.ndarray
	@param z_coordinate: z coordinates as an array
	@type z_coordinate: np.ndarray
	@param bin_number: number of bins
	@type bin_number: int
	@return: Figure
	@rtype: plt.Figure
	"""

	print("Beginning density plot")
	# Compute the radius of each particle from the center (0, 0, 0)
	radius = np.sqrt(x_coordinate ** 2 + y_coordinate ** 2 + z_coordinate ** 2)

	# Define the number of radial bins
	# the particles are distributed in logarithmic fashion

	# calculate the start and end exponent, as the start value in logspace
	# is base ** value_given
	start_exp = math.floor(math.log10(np.min(radius)))
	end_exp = math.ceil(math.log10(np.max(radius)))

	# get logarithmically distributed bins
	bins = np.logspace(start_exp, end_exp, bin_number, base=10.0)
	widths = (bins[1:] - bins[:-1])

	# Compute the histogram of particle counts within radial bins
	hist, bin_edges = np.histogram(radius, bins=bins)

	# Compute the bin radius (middle of each bin), this is a radius on the logarithmic scale
	bin_radius = (bin_edges[:-1] + bin_edges[1:]) / 2

	# Compute the inner and outer radii
	inner_radius = bin_edges[:-1]
	outer_radius = bin_edges[1:]

	# Compute the volume of each bin (spherical shells)
	# QUESTION
	# This stays in the logarithmic scale -> does not influence the particles
	# that are contained in a shell?
	bin_volume = 4. / 3. * np.pi * (outer_radius ** 3 - inner_radius ** 3)

	# Compute the radial number density for each bin
	# take mass of the shell -> count number of partciles * mass
	masses_per_shell = hist * np.min(mass)
	density_per_shell = masses_per_shell / bin_volume

	# Calculate the error of the density in each bin
	# This is the absolute error of each shell
	error_of_density = (np.min(mass) / bin_volume) * np.sqrt(hist)
	#error_of_density = np.divide(np.min(mass), bin_volume, out=np.zeros_like(hist), where=hist != 0, casting="unsafe") * np.sqrt(hist)

	# relative error
	#rel_density_error = np.divide(np.sqrt(hist), hist, out=np.zeros_like(hist), where=hist != 0, casting="unsafe")
	#rel_density_error = np.divide(1, np.sqrt(hist), out=np.zeros_like(hist), where=hist != 0, casting="unsafe")
	#rel_density_error = np.divide(error_of_density, density_per_shell, out=np.zeros_like(error_of_density), where=density_per_shell != 0)

	# fit the data

	# get the total total_mass of the system and give it as an array so it can be used
	tot_mass = np.sum(mass)
	tot_mass_array = np.zeros_like(bin_radius)
	tot_mass_array.fill(tot_mass)

	# set an initial guess for the parameter a, that should be fitted
	initial_guess = 0.1

	# return the fitted parameter a and the covariance
	# maxfev extends the default number of iterations so that the
	# fittted a can be found
	params = curve_fit(hernquist_rho, (bin_radius, tot_mass_array),
						density_per_shell, p0=initial_guess, maxfev=100000)
						#sigma=rel_density_error) # QUESTION how does that work, when implementing it, the value always goes to the inital guess and it divides by 0 as there are empty bins

	# Extract the fitted parameters
	fitted_a = params[0]

	# Generate the fitted curve using the fitted a
	fitted_curve = hernquist_rho((bin_radius, tot_mass), fitted_a)

	# Plot the radial number density profile
	particle_density_plot = plt.figure()
	plt.bar(bin_radius, density_per_shell, width=widths, label="Particle Density")
	plt.errorbar(bin_radius, density_per_shell, xerr=None, yerr=error_of_density, color="orange", ls="none")
	plt.plot(bin_radius, fitted_curve, linestyle='--', color="red", label="Hernquist density")
	plt.title('Radial Number Density Profile')
	plt.legend()
	plt.xlabel('Radius')
	plt.xscale("log")
	plt.ylabel('Number Density')
	plt.yscale("log")

	print("Finished the density plot")

	return particle_density_plot
