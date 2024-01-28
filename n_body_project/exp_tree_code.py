import matplotlib.pyplot as plt
import numpy as np

from n_body_shells import max_radius


# start out with a quadtree to get the idea
# TODO when it is running as it should, expand it to an OctTree
# TODO if possible change it to jitclass (and leave it, when it requires a rewrite of the whole code)

class Point:
	def __init__(self, x_coordinate, y_coordinate, z_coordinate, mass, ident):
		self.x = x_coordinate
		self.y = y_coordinate
		self.z = z_coordinate
		self.mass = mass
		self.ident = ident


		self.x_accel = 0
		self.y_accel = 0
		self.z_accel = 0

	def __str__(self):
		return f"Point {self.ident}: {self.x}, {self.y}, {self.z}"

	def __repr__(self):
		return f"Point ({self.ident} {self.x}, {self.y}, {self.z})"


class OctTreeNode:
	def __init__(self, xmin, ymin, zmin, xmax, ymax, zmax, theta):
		self.xmin = xmin
		self.ymin = ymin
		self.zmin = zmin
		self.xmax = xmax
		self.ymax = ymax
		self.zmax = zmax

		self.xmid = (self.xmin + self.xmax) / 2
		self.ymid = (self.ymin + self.ymax) / 2
		self.zmid = (self.zmin + self.zmax) / 2

		self.total_mass = 0
		self.x_com = 0
		self.y_com = 0
		self.z_com = 0

		self.points = []
		self.children = []

		self.max_items = 1
		self.theta = theta
		self.s = self.xmax - self.xmid


	def add_point(self, point):
		# check if the current node is full
		# if yes, then divide the node
		# self.points is set None in split_node -> it will skip this step
		# when the first point is added and when there are already children
		if self.points is not None and len(self.points) >= self.max_items:
			# if it is, split the points up into different children
			self.split_node()

		# if splitting is not necessary, check if the point can just be
		# appended to the points of the node -> this is only possible
		# if the node has no children.
		# If the node has children, the point has to be added to the correct child
		if self.points is not None: # none marks an internal node, so here the particle is added to an exeternal node
			# TODO calculate the total mass and center of mass of the current node here
			self.points.append(point)
			self.calculate_center_of_mass(point) # from here on, each child node has a center of mass and a total mass
		else:
			self.add_point_to_child(point)

	def split_node(self):

		# this makes the tree recursive
		self.children.append(OctTreeNode(self.xmin, self.ymin, self.zmin, self.xmid, self.ymid, self.zmid, self.theta))
		self.children.append(OctTreeNode(self.xmin, self.ymid, self.zmin, self.xmid, self.ymax, self.zmid, self.theta))
		self.children.append(OctTreeNode(self.xmid, self.ymid, self.zmin, self.xmax, self.ymax, self.zmid, self.theta))
		self.children.append(OctTreeNode(self.xmid, self.ymin, self.zmin, self.xmax, self.ymid, self.zmid, self.theta))
		self.children.append(OctTreeNode(self.xmin, self.ymin, self.zmid, self.xmid, self.ymid, self.zmax, self.theta))
		self.children.append(OctTreeNode(self.xmin, self.ymid, self.zmid, self.xmid, self.ymax, self.zmax, self.theta))
		self.children.append(OctTreeNode(self.xmid, self.ymid, self.zmid, self.xmax, self.ymax, self.zmax, self.theta))
		self.children.append(OctTreeNode(self.xmid, self.ymin, self.zmid, self.xmax, self.ymid, self.zmax, self.theta))

		# divide points onto the nodes
		for current_points in self.points:
			self.add_point_to_child(current_points)

		# reset the point list to 0
		self.points = None # with setting it to None, it makes it an internal node

	def add_point_to_child(self, point):
		# add a given point to the correct child
		# loop over all children, determine to which child it belongs
		# QUESTION this way, the point is only added once, to the node that it fits in first
		for child in self.children:
			if child.xmin <= point.x <= child.xmax and child.ymin <= point.y <= child.ymax and child.zmin <= point.z <= child.zmax:
				child.add_point(point)
				# QUESTION calculate center of mass here for the node
				# after the point has been added to a child,
				# further search is not possible
				break

		# after the children are created, the center of mass is calculated
		self.update_center_of_mass()


	def calculate_center_of_mass(self, point):
		self.total_mass = point.mass
		self.x_com = (point.x * point.mass) / self.total_mass
		self.y_com = (point.y * point.mass) / self.total_mass
		self.z_com = (point.z * point.mass) / self.total_mass

		return self.total_mass, self.x_com, self.y_com, self.z_com

	def update_center_of_mass(self):
		total_mass = 0
		x_com = 0
		y_com = 0
		z_com = 0

		for child in self.children: # FIXME I don't think that works like this for the COM
			total_mass += child.total_mass
			x_com += child.x_com * child.total_mass
			y_com += child.y_com * child.total_mass
			z_com += child.z_com * child.total_mass

		if total_mass > 0:
			# all these values are started out with 0
			# if the total mass does not accumulate to anything of value
			# there is no particle in the node, therefore no center of mass
			# and no total mass of that node can be calculated
			self.total_mass = total_mass
			self.x_com = x_com / total_mass
			self.y_com = y_com / total_mass
			self.z_com = z_com / total_mass

		return self.total_mass, self.x_com, self.y_com

	def calculate_acceleration(self, point):
		"""Calculate the acceleration acting on a particle
		The particle has to be added from the main file
		loop over them to calculate the acceration for each particle"""

		# if the current node is an external node and it is not body b
		# then the node has no children, but contains a point
		if self.points is not None and len(self.points) > 0 and point.ident is not self.points[0].ident:
			rx = point.x - self.points[0].x
			ry = point.y - self.points[0].y
			rz = point.z - self.points[0].z

			abs_r = np.sqrt(rx**2 + ry**2 + rz**2)

			# no multiplication with G, as G = 1

			acc_x = - (self.points[0].mass / abs_r**3) * rx
			acc_y = - (self.points[0].mass / abs_r**3) * ry
			acc_z = - (self.points[0].mass / abs_r**3) * rz

			point.x_accel += acc_x
			point.y_accel += acc_y
			point.z_accel += acc_z

			return point.x_accel, point.y_accel, point.z_accel

		# the point is internal and has children
		# children that do not have children or a particle do not count towards the
		# acceleration
		else:
			# calculate the distance between the particle and nodes center of mass
			dx = point.x - self.x_com
			dy = point.y - self.y_com
			dz = point.z - self.z_com

			d_abs = np.sqrt(dx**2 + dy**2 + dz**2)

			# check if the threshold value is no exceeded by the ratio
			# if d_abs == 0 -> is the its own point
			if d_abs > 0 and self.s / d_abs < self.theta:
				# treat thread like a single body and calculate the acceleration experienced
				# as the acceleration from the nodes COM and total mass
				acc_x = - (self.total_mass / d_abs ** 3) * dx
				acc_y = - (self.total_mass / d_abs ** 3) * dy
				acc_z = - (self.total_mass / d_abs ** 3) * dz

				point.x_accel += acc_x
				point.y_accel += acc_y
				point.z_accel += acc_z

				return point.x_accel, point.y_accel, point.z_accel

			else:
				for child in self.children:
					child.calculate_acceleration(point)


def print_tree(root, level=0):
	"""recursively print the tree structure"""
	if root.points is not None and len(root.points) > 0:
		ide = root.points[0].ident
	else:
		ide = "none"
	print(f"{'  ' * level} level: {level}, mass: {root.total_mass}, number of children: {len(root.children)}, point ID: {ide}")
	for child in root.children:
		print_tree(child, level + 1)


def opening_angles(start, stop, n_angles):
	"""get a list of opening angles for the tree code"""

	if not (0 <= start <= 1 and 0 <= stop <= 1):
		raise ValueError("Start and Stop have to be between 0 and 1")

	angles = np.linspace(start, stop, n_angles)

	angles = list(angles)

	return angles


def tree_code_forces(x_coordinates, y_coordinates, z_coordinates, masses, part_number, theta):
	"""calculate the forces with the tree code and associate them with
	the according radii again"""

	radii, maximal_radius = max_radius(x_coordinates, y_coordinates, z_coordinates)

	# create an empty array for the accelerations
	ax_accel, ay_accel, az_accel = np.zeros_like(x_coordinates), np.zeros_like(y_coordinates), np.zeros_like(z_coordinates)

	print("Start tree generation")

	# QUESTION should other min and max values be used -> this won't give a square otherwise -> influences the subdivision
	octtree = OctTreeNode(xmin=np.min(radii), ymin=np.min(radii), zmin=np.min(radii),
						  xmax=np.max(radii), ymax=np.max(radii), zmax=np.max(radii), theta=theta)

	# get the points into an appropriate format
	point_coordinates = [(x_coordinates[i], y_coordinates[i], z_coordinates[i], masses[i], part_number[i]) for i in range(len(x_coordinates))]

	count = 0
	# built the tree and add the particles
	for point in range(len(point_coordinates)):
		count += 1
		#print(f"Calculating particle {point}")
		point = Point(*point_coordinates[point])
		octtree.add_point(point)

	# calculate the forces for each particle
	for point_cords in range(len(point_coordinates)):
		print(point_cords)
		point = Point(*point_coordinates[point_cords])
		acc = octtree.calculate_acceleration(point)
		ax_accel[point_cords] = point.x_accel
		ay_accel[point_cords] = point.y_accel
		az_accel[point_cords] = point.z_accel

	Fx = ax_accel * masses
	Fy = ay_accel * masses
	Fz = az_accel * masses
	abs_F = np.sqrt(Fx ** 2 + Fy ** 2 + Fz ** 2)

	abs_r = np.sqrt(x_coordinates ** 2 + y_coordinates ** 2 + z_coordinates ** 2)

	accelerations = (ax_accel, ay_accel, az_accel)

	return abs_F, abs_r, accelerations


def repetitive_tree_code(x_coordinates, y_coordinates, z_coordinates, masses, part_number, opening_angles):
	"""calculate the forces with the tree code several times for different opening angles"""

	print("Start Tree Force plotting")

	tree_plot, ax_tree_plot = plt.subplots()

	for theta in opening_angles:
		abs_F, abs_r, accelerations = tree_code_forces(x_coordinates, y_coordinates, z_coordinates, masses, part_number, theta)
		F_abs_sort = np.rec.fromarrays([abs_F, abs_r], dtype=np.dtype([('abs_F', np.float64), ('abs_r', np.float64)]))
		F_abs_sort.sort()

		# extract the sorted arrays to use them for plotting
		r_abs_sort = F_abs_sort.abs_r
		F_abs_sort = F_abs_sort.abs_F

		ax_tree_plot.scatter(r_abs_sort, F_abs_sort, label=f"theta={theta}")

	ax_tree_plot.legend()
	ax_tree_plot.set_xlabel("Absolute Radius")
	ax_tree_plot.set_ylabel("Forces (Tree Code)")
	ax_tree_plot.set_xscale("log")
	ax_tree_plot.set_yscale("log")

	print("End Tree Force plotting")

	return tree_plot





