# libraries

# TODO put failsaves in all function to make sure all particles are processed and in the right order / index

# my functions
from particle_density import *
from n_body_forces import *
from n_body_shells import *
from relaxation_time import *
from exp_tree_code import *
from time_integration import *

# ------------------------------------------------------------------------------

# read in data
# data structure: particle_number, masses, x, y, z, vx, vy, vz, softening, potential
data = np.loadtxt("n_body_project/data/data.txt")

# extract the data
particle_number, mass, x_cord, y_cord, z_cord, vx_vel, vy_vel, vz_vel, softening, potential = data_extraction(data)

# all data in Planck units
# ------------------------------------------------------------------------------
"""
# Task 1: Step 1

# get the density plot
particle_density_plot = particle_density(mass, x_cord, y_cord, z_cord, bin_number=50)
particle_density_plot.savefig("n_body_project/plots/density_profile.png")
# plt.show()
plt.close(particle_density_plot)

# test time integration
#leap_frog(x_cord, y_cord, z_cord, vx_vel, vy_vel, vz_vel, mass, 1, 1)
"""
# ------------------------------------------------------------------------------
# Task 1: Step 2
"""
# calculate half-total_mass radius
hmr, half_mass = half_mass_radius(x_cord, y_cord, z_cord, mass)

# get all the particles within the half total_mass radius
pos = particles_in_hmr(hmr, x_cord, y_cord, z_cord)

# calculate the mean inter particle separation
mean_int = mean_inter_particle_separation(pos)

# get values of different order of magnitude for softening
s = softening_values(mean_int, 3, 3)

# calculate the forces brute force and plot their dependence on the softening
soft_plot, rad_plot = softening_plot(s, particle_number, mean_int, mass, x_cord, y_cord, z_cord, softening)
soft_plot.savefig("n_body_project/plots/softening.png")
#fig, ax = rad_plot
#fig.savefig("n_body_project/plots/rad.png")
# plt.show()
plt.close(soft_plot)
#plt.close(fig)

# calculate the forces depending on the shell
shell_ind, shell_bound, shell_thickness = split_data_into_shells(x_cord, y_cord, z_cord, 1000)
rhoos = rho_of_shell(shell_ind, shell_bound, mass)
masses_of_shells = shell_masses(rhoos, shell_bound, shell_thickness)

forces_on_particles, radii = forces_when_looking_at_shells(masses_of_shells, shell_bound, x_cord, y_cord, z_cord, mass)

# sort the forces and the radii
F_abs_sort = np.rec.fromarrays([forces_on_particles, radii], dtype=np.dtype([('F_abs', np.float64), ('r_abs', np.float64)]))  # find the permutation to sort r_abs fro lowest to highest
F_abs_sort.sort()

r_abs_sort = F_abs_sort.r_abs
F_abs_sort = F_abs_sort.F_abs

fig, ax = rad_plot

ax.scatter(radii, forces_on_particles, color="green", edgecolors='black', label=f'Shell forces')
ax.legend()
fig.savefig("n_body_project/plots/rad.png")
plt.close(fig)

# QUESTION we take the half mass radius for this, should we also take just the number of
# particles included in the half mass radius or all of them
rel_time, cross_time = relaxation_time(len(particle_number), half_mass, hmr)
print(f"relaxation time: {rel_time}, crossing time: {cross_time}")

"""
#-------------------------------------------------------------------------------
# Task 2: Tree code

octtree = OctTreeNode(xmin=np.min(x_cord), ymin=np.min(y_cord), zmin=np.min(z_cord), xmax=np.max(x_cord), ymax=np.max(y_cord), zmax=np.max(z_cord))
point_coordinates = [(x_cord[i], y_cord[i], z_cord[i], mass[i], particle_number[i]) for i in range(len(x_cord))]
point_coordinates = point_coordinates[::100]
print(len(point_coordinates))
print(point_coordinates)

count = 0
for point in range(len(point_coordinates)):
	count += 1
	print(f"Calculating particle {point}")
	point = Point(*point_coordinates[point])
	octtree.add_point(point)

print(count)
print(octtree.total_mass)






