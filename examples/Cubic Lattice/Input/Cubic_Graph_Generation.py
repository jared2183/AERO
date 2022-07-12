from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#Start with a single cell
h = 7.0
w = 7.0
l = 7.0
H = 3
W = 3
L = 3
vertices = []
origin = np.array([0,0,0])

#Generate vertices
# 1. 5 on bottom in X shape
# 2. 4 in middle in diamond shape
# 3. 5 on top in X shape

#Bottom layer
for k in range(H+1):
	for j in range(L+1):
		for i in range(W+1):
			vertices.append(origin)
			origin = origin + [w,0,0]
		origin[0] = 0
		#Change origin
		origin = origin + [0,l,0]
	origin[0] = origin[1] = 0
	#Change origin
	origin = origin +  [0,0,h]

vertices = np.unique(vertices,axis=0)
#import ipdb; ipdb.set_trace()
edges = []


for index_master,vertex_master in enumerate(vertices):
	for index_slave,vertex_slave in enumerate(vertices):
		if index_slave != index_master:
			if abs(np.linalg.norm(vertex_master-vertex_slave)) == h:
				edges.append([index_master,index_slave])
				plt.plot([vertex_master[0],vertex_slave[0]],[vertex_master[1],vertex_slave[1]],[vertex_master[2],vertex_slave[2]])

for vertex in vertices:
	ax.scatter(*vertex)
plt.show()

edges = np.array(edges)
edges.sort(axis=1)
edges = np.unique(edges,axis=0)

np.savetxt('cubic_edges.csv',edges,fmt='%i',delimiter=',')
np.savetxt('cubic_vertices.csv',vertices,fmt='%10.6f',delimiter=',')


