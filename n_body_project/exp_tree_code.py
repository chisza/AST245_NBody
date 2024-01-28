import matplotlib.pyplot as plt
import numpy as np


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

	def __str__(self):
		return f"Point {self.ident}: {self.x}, {self.y}, {self.z}"

	def __repr__(self):
		return f"Point ({self.ident} {self.x}, {self.y}, {self.z})"


class OctTreeNode:
	def __init__(self, xmin, ymin, zmin, xmax, ymax, zmax):
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
		self.theta = 0.5

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
		self.children.append(OctTreeNode(self.xmin, self.ymin, self.zmin, self.xmid, self.ymid, self.zmid))
		self.children.append(OctTreeNode(self.xmin, self.ymid, self.zmin, self.xmid, self.ymax, self.zmid))
		self.children.append(OctTreeNode(self.xmid, self.ymid, self.zmin, self.xmax, self.ymax, self.zmid))
		self.children.append(OctTreeNode(self.xmid, self.ymin, self.zmin, self.xmax, self.ymid, self.zmid))
		self.children.append(OctTreeNode(self.xmin, self.ymin, self.zmid, self.xmid, self.ymid, self.zmax))
		self.children.append(OctTreeNode(self.xmin, self.ymid, self.zmid, self.xmid, self.ymax, self.zmax))
		self.children.append(OctTreeNode(self.xmid, self.ymid, self.zmid, self.xmax, self.ymax, self.zmax))
		self.children.append(OctTreeNode(self.xmid, self.ymin, self.zmid, self.xmax, self.ymid, self.zmax))

		# divide points onto the nodes
		for current_points in self.points:
			self.add_point_to_child(current_points)

		# reset the point list to 0
		self.points = None # with setting it to None, it makes it an internal node

		# calculate the center of mass for each node (child) -> each child has a center of mass and a total mass
		#self.update_center_of_mass()

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

		for child in self.children:
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


# Implementation
def plot_quadtree(node, ax=None):
	if ax is None:
		_, ax = plt.subplots()
#
	ax.plot([node.xmin, node.xmin, node.xmax, node.xmax, node.xmin],
			[node.ymin, node.ymax, node.ymax, node.ymin, node.ymin], color='black')
#
	for child in node.children:
		plot_quadtree(child, ax)

	if node.points is not None:
		for point in node.points:
			ax.plot(point.x, point.y, 'ro')
#
	ax.set_xlim(0, 100)
	ax.set_ylim(0, 100)
	ax.set_aspect('equal', adjustable='box')
#
#
quadtree = OctTreeNode(xmin=0, ymin=0, xmax=100, ymax=100, zmin=0, zmax=100)
points_list = [(10, 10, 3, 1, 1), (20, 20, 4, 1, 1), (5, 5, 5, 1, 1), (90, 90, 6, 1, 1), (15, 15, 7, 1, 1), (30, 30, 8, 1, 1),
			   (1, 1, 1, 1, 1), (1, 40, 1, 1, 1), (20, 40, 1, 1, 1), (1, 49, 1, 1, 1)]
#
for point_coords in points_list:
	point = Point(*point_coords)
	quadtree.add_point(point)

print(quadtree.total_mass)
#
plot_quadtree(quadtree)
plt.show()
#



