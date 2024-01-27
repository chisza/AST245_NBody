import matplotlib.pyplot as plt
import numpy as np


# start out with a quadtree to get the idea
# TODO when it is running as it should, expand it to an OctTree
# TODO if possible change it to jitclass (and leave it, when it requires a rewrite of the whole code)

class Point:
	def __init__(self, x_coordinate, y_coordinate, mass):
		self.x = x_coordinate
		self.y = y_coordinate
		self.mass = mass
		#self.id = id

		self.x_accel = 0
		self.y_accel = 0

	def __str__(self):
		return f"Point {self.id}: {self.x}, {self.y}"

	def __repr__(self):
		return f"Point ({self.id} {self.x}, {self.y})"


class QuadTreeNode:
	def __init__(self, xmin, ymin, xmax, ymax):
		self.xmin = xmin
		self.ymin = ymin
		self.xmax = xmax
		self.ymax = ymax

		self.total_mass = 0
		self.x_com = 0
		self.y_com = 0

		self.points = []
		self.children = []
		self.max_items = 1
		self.theta = 0.5

	def __str__(self, level=0):
		ret = " " * level + str(self.x_com) + '\n'
		for child in self.children:
			ret += child.__str__(level + 1)

		return ret


	def add_point(self, point):
		# check if the current node is full
		# if yes, then divide the node
		if self.points is not None and len(self.points) >= self.max_items:
			# if it is, split the points up into different children
			self.split_node()

		# if splitting is not necessary, check if the point can just be
		# appended to the points of the node -> this is only possible
		# if the node has no children.
		# If the node has children, the point has to be added to the correct child
		if self.children:
			self.add_point_to_child(point)
			print("we have children")
		else:
			self.points.append(point)
			self.calculate_center_of_mass(point)

	def calculate_center_of_mass(self, point):
		# just add the new point coordinates and mass to the existing center of mass coordinates and mass
		total_mass = self.total_mass + point.mass
		x_com = (self.x_com * self.total_mass + point.x * point.mass) / total_mass
		y_com = (self.y_com * self.total_mass + point.y * point.mass) / total_mass

		self.total_mass = total_mass
		self.x_com = x_com
		self.y_com = y_com

	def split_node(self):
		# calculate the midpoints for the child nodes
		child_xmid = (self.xmin + self.xmax) / 2
		child_ymid = (self.ymin + self.ymax) / 2

		# this makes the tree recursive
		self.children.append(QuadTreeNode(self.xmin, self.ymin, child_xmid, child_ymid))
		self.children.append(QuadTreeNode(self.xmin, child_ymid, child_xmid, self.ymax))
		self.children.append(QuadTreeNode(child_xmid, child_ymid, self.xmax, self.ymax))
		self.children.append(QuadTreeNode(child_xmid, self.ymin, self.xmax, child_ymid))

		# divide points onto the nodes
		for current_points in self.points:
			self.add_point_to_child(current_points)

		# reset the point list to 0
		self.points = []

		# after splitting the node, each child needs the center of total_mass coordinates
		# and the total total_mass calculated
		for child in self.children:
			pass
			#child.calculate_center_of_mass()

	def update_center_of_mass(self):
		# check if the node has children
		if self.children:
			total_mass = 0
			x_com = 0
			y_com = 0

			for child in self.children:
				child.update_center_of_mass()
				total_mass += child.total_mass
				x_com += child.x_com * child.total_mass
				y_com += child.y_com * child.total_mass

			if total_mass > 0:
				self.total_mass = total_mass
				self.x_com = x_com / total_mass
				self.y_com = y_com / total_mass



	def add_point_to_child(self, point):
		# add a given point to the correct child
		# loop over all children, determine to which child it belongs
		for child in self.children:
			if child.xmin <= point.x <= child.xmax and child.ymin <= point.y <= child.ymax:
				child.add_point(point)
				# QUESTION calculate center of mass here for the node
				# after the point has been added to a child,
				# further search is not possible
				break

	# determine if a node is far enough from the root particle that it can be approximated
	# with the center of mass for the force calculation
	# do this outside, go through the tree and find each point, access its node
	# and go from there
	def force_calculation(self, point):
		vector_radius = np.array((point.x - self.x_com, point.y - self.y_com))
		print(vector_radius)
		abs_radius = np.linalg.norm(vector_radius)
		print(abs_radius)
		print(f"abs_radius type{type(abs_radius)}")
		# check distance
		# FIXME This is not working
		if (self.xmax - self.xmin) / abs_radius < self.theta:
			# if distance okay
			# force calculation with com
			print("This node is far away enough")
			# else:
			# go recursivley through children
			if self.children:
				print("Hello")

		else:
			if self.points is not None and len(self.points) == 1:
				# when there are no points in the node, we do not need to calculate it
				# as the node is empty -> no force
				# calculate the distance between the root particle
				# and the current particle
				G = 1.
				#vector_radius = np.array((point.x - self.x_com, point.y - self.y_com))
				#abs_radius = np.sqrt(vector_radius ** 2)
				x_acceleration = - (G * self.total_mass) / abs_radius ** 2 * vector_radius[0]
				y_acceleration = - (G * self.total_mass) / abs_radius ** 2 * vector_radius[1]
				print("calc accel")
				point.x_accel += x_acceleration
				point.y_accel += y_acceleration


# Implementation
# def plot_quadtree(node, ax=None):
# 	if ax is None:
# 		_, ax = plt.subplots()
# #
# 	ax.plot([node.xmin, node.xmin, node.xmax, node.xmax, node.xmin],
# 			[node.ymin, node.ymax, node.ymax, node.ymin, node.ymin], color='black')
# #
# 	for child in node.children:
# 		plot_quadtree(child, ax)
# #
# 	for point in node.points:
# 		ax.plot(point.x, point.y, 'ro')
# #
# 	ax.set_xlim(0, 100)
# 	ax.set_ylim(0, 100)
# 	ax.set_aspect('equal', adjustable='box')
# #
# #
# quadtree = QuadTreeNode(xmin=0, ymin=0, xmax=100, ymax=100)
# points_list = [(10, 10, 3), (20, 20, 4), (5, 5, 5), (90, 90, 6), (15, 15, 7), (30, 30, 8), (1, 1, 1), (1, 40, 1), (20, 40, 1), (1, 49, 1)]
# #
# for point_coords in points_list:
# 	point = Point(*point_coords)
# 	quadtree.add_point(point)
# #
# plot_quadtree(quadtree)
# plt.show()
#



