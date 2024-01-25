# libraries

# TODO put failsaves in all function to make sure all particles are processed and in the right order / index

# my functions
from particle_density import *
from n_body_forces import *
from n_body_shells import *
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

# ------------------------------------------------------------------------------
# Task 1: Step 2

# calculate half-total_mass radius
hmr = half_mass_radius(x_cord, y_cord, z_cord, mass)

# get all the particles within the half total_mass radius
pos = particles_in_hmr(hmr, x_cord, y_cord, z_cord)

# calculate the mean inter particle separation
mean_int = mean_inter_particle_separation(pos)

# get values of different order of magnitude for softening
s = softening_values(mean_int, 1, 1)

# calculate the forces brute force and plot their dependence on the softening
soft_plot, rad_plot = softening_plot(s, particle_number, mean_int, mass, x_cord, y_cord, z_cord, softening)
soft_plot.savefig("n_body_project/plots/softening.png")
#rad_plot.savefig("n_body_project/plots/rad.png")
# plt.show()
plt.close(soft_plot)
#plt.close(rad_plot)

# calculate the forces depending on the shell
shell_ind, shell_bound, shell_thickness = split_data_into_shells(x_cord, y_cord, z_cord, 50)
rhoos = rho_of_shell(shell_ind, shell_bound, mass, particle_number)
masses_of_shells = shell_masses(rhoos, shell_bound, shell_thickness)

forces_on_particles, radii = forces_when_looking_at_shells(masses_of_shells, shell_bound, x_cord, y_cord, z_cord)

# sort the forces and the radii
F_abs_sort = np.rec.fromarrays([forces_on_particles, radii], dtype=np.dtype([('F_abs', np.float32), ('r_abs', np.float32)]))  # find the permutation to sort r_abs fro lowest to highest
F_abs_sort.sort()

r_abs_sort = F_abs_sort.r_abs
F_abs_sort = F_abs_sort.F_abs

fig, ax = rad_plot

ax.scatter(radii[::100], forces_on_particles[::100], color="green", edgecolors='black', label=f'Shell forces')
ax.legend()
fig.savefig("n_body_project/plots/rad.png")
plt.close(fig)

"""

#-------------------------------------------------------------------------------
# Task 2: Tree code

def plot_quadtree(node, ax=None):
	if ax is None:
		_, ax = plt.subplots()

	ax.plot([node.xmin, node.xmin, node.xmax, node.xmax, node.xmin],
			[node.ymin, node.ymax, node.ymax, node.ymin, node.ymin], color='black')

	for child in node.children:
		plot_quadtree(child, ax)

	for point in node.points:
		ax.plot(point.x, point.y, 'ro')

	ax.set_xlim(0, np.max(x_cord))
	ax.set_ylim(0, np.max(y_cord))
	ax.set_aspect('equal', adjustable='box')


quadtree = QuadTreeNode(xmin=0, ymin=0, xmax=np.max(x_cord), ymax=np.max(y_cord))
point_coordinates = [(x_cord[i], y_cord[i], mass[i]) for i in range(len(x_cord))]
#point_coordinates = point_coordinates[::10]
print(len(point_coordinates))
print(point_coordinates)

count = 0
for point in range(len(point_coordinates)):
	count += 1
	print(f"Calculating particle {point}")
	point = Point(*point_coordinates[point])
	quadtree.add_point(point)

print(count)

plot_quadtree(quadtree)
plt.show()






