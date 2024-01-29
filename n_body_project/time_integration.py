import numpy as np
from matplotlib import pyplot as plt
from numba import jit, objmode

from n_body_forces import *

from n_body_forces import *
from relaxation_time import *


@jit(nopython=True)
def leap_frog(x_cord, y_cord, z_cord, vel_x, vel_y, vel_z, masses, softening, h):
	"""Using the leap frog method to calculate the progression of the system over time

	@param r0: the initial position as array (x, y, z)
	@param v0: the initial velocity as arry (x, y, z)
	@param masses: the masses of the particles
	@param softening: the softening for the calculation
	@param h: the timestep
	"""
	print("here")
	# drift motion: calculate the half position
	r_step_half_x = x_cord + 0.5 * h * vel_x
	r_step_half_y = y_cord + 0.5 * h * vel_y
	r_step_half_z = z_cord + 0.5 * h * vel_z

	print("here2")

	# calculate all the accelerations at once
	acceleration = direct_force_calculation(masses, r_step_half_x, r_step_half_y, r_step_half_z, softening)[2]

	# kick motion: calculate the new velocity
	v_step_one_x = vel_x + h * acceleration[0]
	v_step_one_y = vel_y + h * acceleration[1]
	v_step_one_z = vel_z + h * acceleration[2]

	# full drift motion: calculate the new position
	r_step_one_x = r_step_half_x + 0.5 * h * v_step_one_x
	r_step_one_y = r_step_half_y + 0.5 * h * v_step_one_y
	r_step_one_z = r_step_half_z + 0.5 * h * v_step_one_z

	print("this runs, if it is correct, is another question")

	new_position = (r_step_one_x, r_step_one_y, r_step_one_z)
	new_velocity = (v_step_one_x, v_step_one_y, v_step_one_z)

	return new_position, new_velocity

def repetitive_leapfog(x_cord, y_cord, z_cord, vel_x, vel_y, vel_z, masses, softening, time, h):
	"""Repeat the leap frog algorithm several times

	@param time: the total time over which the integration should run
	@param h: Stepsize choosen for the integration

	"""

	stuff_for_plotting = []

	steps = math.floor(time / h)

	for i in range(steps):

		# calculate the new position and the new velocity
		new_position, new_velocity = leap_frog(x_cord, y_cord, z_cord, vel_x, vel_y, vel_z, masses, softening, h)

		# assign the new position and the new velocity
		x_cord = new_position[0]
		y_cord = new_position[1]
		z_cord = new_position[2]

		vel_x = new_velocity[0]
		vel_y = new_velocity[1]
		vel_z = new_velocity[2]

		# get a list of all the values for plotting
		stuff_for_plotting.append((x_cord, y_cord, z_cord))

		print(f"done {i}")

		# TODO return something for plotting
	return x_cord, y_cord, z_cord


def relaxation_magnitude(x_cord, y_cord, z_cord, vel_x, vel_y, vel_z, masses, softening, rep, t_cross):
	"""Let leapfrog run for several timesteps and see how the time scale changes"""

	time = rep * t_cross
	time = math.ceil(time)

	h = t_cross

	steps = time / h
	steps = math.floor(steps)

	for i in range(0, time, steps):

		for j in range(len(softening)):

			soft = int(j)

			new_position, new_velocity = leap_frog(x_cord, y_cord, z_cord, vel_x, vel_y, vel_z, masses, soft, h)

			# assign the new position and the new velocity
			x_cord = new_position[0]
			y_cord = new_position[1]
			z_cord = new_position[2]

			vel_x = new_velocity[0]
			vel_y = new_velocity[1]
			vel_z = new_velocity[2]

			# calculate half mass radius
			hmr, half_mass = half_mass_radius(x_cord, y_cord, z_cord, masses)

			# particles in the half masss radius
			part_in_hmr = particles_in_hmr(hmr, x_cord, y_cord, z_cord)

			# calculate relaxation timescale
			relaxation_time_system, crossing_time_system = relaxation_time(part_in_hmr, half_mass, hmr)

			print(f"iteration: {i}, softening: {softening[j]}, relaxation time: {relaxation_time_system}, crossing time: {crossing_time_system}, half mass radius: {half_mass}")







