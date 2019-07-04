# Python program print Eulerian Trail in a given Eulerian or Semi-Eulerian Graph

from collections import defaultdict, OrderedDict, Counter
import copy
import numpy as np
import matplotlib
from pathos.multiprocessing import ProcessingPool as Pool
import sys

if sys.platform == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")

plot_results = True

#This class represents an undirected graph using adjacency list representation
class Graph:
    def __init__(self,vertices,vertice_vals):
        self.V= vertices #No. of vertices
        self.E = []
        self.graph = defaultdict(list) # default dictionary to store graph
        self.intersections = defaultdict(list)
        self.finished= False
        self.vert_vals = vertice_vals
        self.color = None
        self.edge_count =0
        self.colinear_vals = []
        self.kill_list = []
        self.bad_angle_count = 0

    # This function adds edge u-v to graph
    def addEdge(self,u,v):
        if v not in self.graph[u]:
            #Only add edge if it satisfies angle requirement
            print_line = self.vert_vals[v]-self.vert_vals[u]
            print_line = self.vert_vals[1]-self.vert_vals[6]
            mag_line = np.sqrt(print_line[0]**2+print_line[1]**2+print_line[2]**2)
            vertical_line = np.array([0,0,1])
            dot_p = np.sum(vertical_line*print_line)
            print_angle = np.arccos(dot_p/mag_line)/np.pi*180

            if print_angle < 150:
                self.edge_count += 1
                self.graph[u].append(v)
                self.graph[v].append(u)
            else:
                self.graph[v].append(u)

    # This function removes edge u-v from graph    
    def rmvEdge(self, u, v):
        edge_removed = False
        for index, key in enumerate(self.graph[u]):
            if key == v:
                self.graph[u].pop(index)
                edge_removed = True
                self.edge_count -=1
        for index, key in enumerate(self.graph[v]):
            if key == u:
                self.graph[v].pop(index)
                self.edge_count -=1
        return edge_removed
 
    # A DFS based function to count reachable vertices from v
    def DFSCount(self, v, visited):
        count = 1
        visited[v] = True
        for i in self.graph[v]:
            if visited[i] == False:
                count = count + self.DFSCount(i, visited)         
        return count
 
    def isValidAngle(self, u, v):
        # Check if moving from u to v is a valid printing move for embedded in terms of angle
        # from: https://math.stackexchange.com/questions/974178/how-to-calculate-the-angle-between-2-vectors-in-3d-space-given-a-preset-function
        
        print_line = self.vert_vals[v]-self.vert_vals[u]
        mag_line = np.sqrt(print_line[0]**2+print_line[1]**2+print_line[2]**2)
        vertical_line = np.array([0,0,1])
        dot_p = np.sum(vertical_line*print_line)
        print_angle = np.arccos(dot_p/mag_line)/np.pi*180
        
        if print_angle > 150:
            #import matplotlib.pyplot as plt
            #from mpl_toolkits.mplot3d import Axes3D
            #fig = plt.figure()
            #ax = fig.add_subplot(111, projection='3d')
            #ax.plot([0,vertical_line[0]],[0,vertical_line[1]],[0,vertical_line[2]])
            #ax.plot([0,print_line[0]],[0,print_line[1]],[0,print_line[2]])
            #ax.set_aspect('equal')
            #ax.set_xlim([-6,6])
            #ax.set_ylim([-6,6])
            #ax.set_zlim([-6,6])
            #plt.show()
            #import ipdb; ipdb.set_trace()
            self.bad_angle_count +=1
            return False
        
        else:
            return True

    # The function to check if edge u-v can be considered as next edge in Euler Tour
    def isValidNextEdge(self, u, v):
        # If the edge u-v is the correct angle, it is valid in one of the following two cases:
        # 1) If v is the only adjacent vertex of u
        if len(self.graph[u]) == 1:
            return True
        else:
            '''
             2) If there are multiple adjacents, then u-v is not a bridge
                 Do following steps to check if u-v is a bridge
            '''
            # 2.a) count of vertices reachable from u  
            visited =[False]*(self.V)
            count1 = self.DFSCount(u, visited)
 
            # 2.b) Remove edge (u, v) and after removing the edge, count vertices reachable from u
            edge_removed = self.rmvEdge(u, v)
            visited =[False]*(self.V)
            count2 = self.DFSCount(u, visited)
 
            # 2.c) Add the edge back to the graph
            if edge_removed:
                self.addEdge(u,v)
 
            # 2.d) If count1 is greater, then edge (u, v) is a bridge
            if count1 > count2 and not self.finished:
                self.rmvEdge(u, v)
                self.path_order.append(v)
                self.printEulerUtil(v)

            return False if count1 > count2 else True
 
    # Print Euler tour starting from vertex u
    def printEulerUtil(self, u):
        #Recur for all the vertices adjacent to this vertex
        iter_graph = copy.copy(self.graph[u])
        for v in iter_graph:
            #If edge u-v is not removed and it's a a valid next edge
            if self.isValidNextEdge(u, v) and not self.finished:
                self.rmvEdge(u, v)
                self.path_order.append(v)
                self.printEulerUtil(v)
        self.finished = True

 
    '''The main function that print Eulerian Trail. It first finds an odd
   degree vertex (if there is any) and then calls printEulerUtil()
   to print the path '''
    def printEulerTour(self,u=0):
        # u is the starting vertex
        self.path_order = [u]
        self.printEulerUtil(u)
        self.finished = False
        return self.path_order
    
    def addTriangle(self, F):
        ### EXPERIMENTAL ###
        #Add all edges of a triangle object iff it doesn't block access to another edge
        for t in F:
            if not [t[0],t[1]] in self.E or not [t[1],t[0]] in self.E:
                self.E.append([t[0],t[1]])
            if not [t[0],t[2]] in self.E or not [t[1],t[2]] in self.E:
                self.E.append([t[0],t[2]])
            if not [t[1],t[2]] in self.E or not [t[1],t[2]] in self.E:
                self.E.append([t[1],t[2]])

    def lineCheck(self,line1, line2, recursive_stop = False, recursive_stop_2 = True):
        #Check if printing line1 will prevent printing line2
        def onSegment(p,q,r):
            #Check if point r is on line segment pq
            def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
            #Function for comparing floats
            #https://stackoverflow.com/questions/5595425/what-is-the-best-way-to-compare-floats-for-almost-equality-in-python
                return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
            
            if not isclose((q[0]-p[0]),0):
                line_slope = (q[1]-p[1])/(q[0]-p[0])
                constant = p[1]-line_slope*p[0]
                
                if isclose(r[1]-line_slope*r[0],constant):
                    #Point is on line
                    #Check if point is in segment
                    if ((r[0] <= max(p[0], q[0]) or isclose(r[0],max(p[0], q[0]))) and (r[0] >= min(p[0], q[0]) or isclose(r[0],min(p[0], q[0]))) and
                        (r[1] <= max(p[1], q[1]) or isclose(r[1],max(p[1], q[1]))) and (r[1] >= min(p[1], q[1]) or isclose(r[1],min(p[1], q[1])))):
                        return True
                    else:
                        return False
                else:
                    return False
            else:
                #Line being compared is vertical (infinite slope)
                #Check if point is in segment
                if isclose(r[0],p[0]) and (r[1] <= max(p[1], q[1]) or isclose(r[1],max(p[1], q[1]))) and (r[1] >= min(p[1], q[1]) or isclose(r[1],min(p[1], q[1]))):
                    return True
                else:
                    return False

        xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
        ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

        def det(a, b):
            return a[0] * b[1] - a[1] * b[0]

        div = det(xdiff, ydiff)

        def py_ang(line1, line2):
            """ Returns the angle in radians between vectors 'v1' and 'v2' """
            v1 = [line1[1][0]-line1[0][0],line1[1][1]-line1[0][1]]
            v2 = [line2[1][0]-line2[0][0],line2[1][1]-line2[0][1]]
            cosang = np.dot(v1, v2)
            sinang = np.linalg.norm(np.cross(v1, v2))
            return np.arctan2(sinang, cosang)

        def isclose(a, b, rel_tol=1e-09, abs_tol=1e-09):
            #Function for comparing floats
            #https://stackoverflow.com/questions/5595425/what-is-the-best-way-to-compare-floats-for-almost-equality-in-python
            #return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
            return abs(a-b) <= abs_tol

        ### PARALLEL CASE ###
        if isclose(div,0): #Lines are parallel
            #Check if line1 overlaps line2, including endpoints
            # if 0 endpoints are on segment they must not overlap
            # If one endpoint is on segment they must be opposing
            # if two endpoints are on segment they must be overlapping
            if (onSegment(line2[0],line2[1],line1[0]) and onSegment(line2[0],line2[1],line1[1])) or (onSegment(line1[0],line1[1],line2[0]) and onSegment(line1[0],line1[1],line2[1])):
                # Line1 overlaps line2
                # Solution assumes no intersection other plane when overlapping in XY plane
                # GOOD TO PRINT LINE1 IF LOWER IN Z THAN LINE2 (LINE 2 STILL ACCESSIBLE)
                # Since they are parallel as viewed from XY plane (in same plane), whichever average Z is more will be higher
                average_line1 = (line1[0][2] + line1[1][2])/2.0
                average_line2 = (line2[0][2] + line2[1][2])/2.0
                if average_line2 > average_line1:
                    #print 'lines are parallel and overlap, but line2 still ACCESSIBLE'
                    # GOOD to print
                    return True
                else:
                    #print 'lines are parallel and overlap and line2 is NOT ACCESSIBLE'
                    # NOT GOOD to print
                    return False
            elif onSegment(line2[0],line2[1],line1[0]) or onSegment(line2[0],line2[1],line1[1]):
                #print "lines are parallel but facing opposing directions"
                # GOOD to print
                return True
            else:
                # Lines are parallel but don't overlap
                # PRINTING LINE1 WILL NOT STOP LINE2 FROM BEING PRINTED
                #print 'lines are parallel and dont overlap'
                # GOOD to print
                return True

        else: # Lines arent parallel, thus they intersect
            d = (det(*line1), det(*line2))
            x = det(d, xdiff) / div
            y = det(d, ydiff) / div
            intersect_point = [x,y]
            # Do the segments intersect if intersection is on either segment
            if onSegment(line1[0],line1[1],intersect_point) and onSegment(line2[0],line2[1],intersect_point):
                # Segments intersect
                # First Check if they share the same node (in XY)
                all_points = np.array([line1[0][:2],line1[1][:2],line2[0][:2],line2[1][:2]])
                num_points = len(all_points)
                unique_points,unique_indices = np.unique(all_points.round(decimals=6),axis=0,return_inverse=True)
                num_unique_points = len(unique_points)
                
                if num_unique_points < num_points: #There must have been duplicates
                    all_points_inclZ = np.array([line1[0],line1[1],line2[0],line2[1]])
                    unique_points_inclZ = np.unique(all_points_inclZ.round(decimals=6),axis=0) 

                    #Points share a node in 3D
                    if len(unique_points_inclZ) < len(all_points_inclZ):
                        #Used to calculate the locations of duplicates
                        temp_uniques , temp_indices, temp_counts = np.unique(unique_indices,return_index=True,return_counts=True)
                        duplicate_indices = np.argwhere(unique_indices == unique_indices[temp_indices[int(np.argwhere(temp_counts>1))]])
                        non_duplicate_indices = [item for item in range(4) if item not in duplicate_indices]
                        all_points_z = np.array([line1[0][2],line1[1][2],line2[0][2],line2[1][2]])

                        # Essentially if the first line is lower at the intersecting node, then it is good to print
                        if float(all_points_z[non_duplicate_indices[0]]) <= float(all_points_z[non_duplicate_indices[1]]):
                            return True
                        else:
                            # Calculate angle between the lines, if too close, print lower one first
                            # Arrange lines to start with duplciate indice before comparing angle
                            line1_test = [all_points[int(duplicate_indices[0])],all_points[int(non_duplicate_indices[0])]]
                            line2_test = [all_points[int(duplicate_indices[1])],all_points[int(non_duplicate_indices[1])]]
                            if abs(py_ang(line1_test, line2_test)) > np.pi/6.0:
                                return True
                            else:
                                return False

                    else: 
                        #Used to calculate the locations of duplicates
                        temp_uniques , temp_indices, temp_counts = np.unique(unique_indices,return_index=True,return_counts=True)
                        duplicate_indices = np.argwhere(unique_indices == unique_indices[temp_indices[int(np.argwhere(temp_counts>1))]])

                        all_points_z = np.array([line1[0][2],line1[1][2],line2[0][2],line2[1][2]])

                        # Essentially if the first line is lower at the intersecting node, then it is good to print
                        if float(all_points_z[duplicate_indices[0]]) <= float(all_points_z[duplicate_indices[1]]):
                            return True
                        else:
                            return False
                
                else:
                    # GOOD TO PRINT LINE1 IF LOWER IN Z THAN LINE2 (LINE 2 STILL ACCESSIBLE)
                    # Lines intersect and do not share a node
                    # Cannot just use average, since not in same plane
                    # Line2 is higher than line1 if Z value is higher at XY coordinate of intersection
                    try:    
                        if not isclose(line1[1][0]-line1[0][0],0):
                            fraction_of_line_at_intersect_line1 = abs(intersect_point[0]-line1[0][0])/abs(line1[1][0]-line1[0][0])
                        else:
                            fraction_of_line_at_intersect_line1 = abs(intersect_point[1]-line1[0][1])/abs(line1[1][1]-line1[0][1])
                        if not isclose(line2[1][0]-line2[0][0],0):
                            fraction_of_line_at_intersect_line2 = abs(intersect_point[0]-line2[0][0])/abs(line2[1][0]-line2[0][0])
                        else:    
                            fraction_of_line_at_intersect_line2 = abs(intersect_point[1]-line2[0][1])/abs(line2[1][1]-line2[0][1])
                    except:
                        import ipdb;ipdb.set_trace()
                    height_at_intersect_line1 = fraction_of_line_at_intersect_line1*(line1[1][2]-line1[0][2])+line1[0][2]
                    height_at_intersect_line2 = fraction_of_line_at_intersect_line2*(line2[1][2]-line2[0][2])+line2[0][2]

                   
                    if height_at_intersect_line2 > height_at_intersect_line1:
                        #Line2 passes over line1, thus it is fine to print line1
                        # GOOD to print
                        return True
                    elif isclose(height_at_intersect_line1,height_at_intersect_line2):
                        #Lines either share the same node or could have collied in mid-air
                        #Either way it should be fine to print, no benefit either way
                        # GOOD to print
                        return True
                    else:
                        #Line2 passes under line1, thus printing line1 would block access to line2
                        # NOT GOOD to print
                        return False
            else:
                if not recursive_stop:
                    # Segments don't intersect
                    #print "segments dont intersect"
                    # should add section to see if lines are close to intersecting though
                    nozzle_radius = 0.72/2

                    #Apply translation to first line perpendicular in both directions and check for collisions
                    if not isclose(line1[1][0]-line1[0][0],0) and not isclose(line1[1][1]-line1[0][1],0):
                        line1_slope = (line1[1][1]-line1[0][1])/(line1[1][0]-line1[0][0])
                        line1_perpendicular = -1/line1_slope
                        x_move = np.sqrt(1/(1+line1_perpendicular**2))
                        y_move = x_move*line1_perpendicular
                
                    elif isclose(line1[1][0]-line1[0][0],0):
                        x_move = 1
                        y_move = 0

                    else:
                        x_move = 0
                        y_move = 1

                    line1_temp = np.array(line1)

                    # Apply first move:
                    line1_temp[:,0] += x_move*nozzle_radius
                    line1_temp[:,1] += y_move*nozzle_radius

                    check_1 = self.lineCheck(line1_temp, line2, recursive_stop = True)
                    #if not self.lineCheck(line1_temp, line2, recursive_stop = True):
                    #return False

                    # Apply second move:
                    line1_temp[:,0] -= 2*x_move*nozzle_radius
                    line1_temp[:,1] -= 2*y_move*nozzle_radius

                    check_2 = self.lineCheck(line1_temp, line2, recursive_stop = True)
                    #if not self.lineCheck(line1_temp, line2, recursive_stop = True):
                    #return False

                    #Apply translation to first line perpendicular in both directions and check for collisions
                    if not isclose(line2[1][0]-line2[0][0],0) and not isclose(line2[1][1]-line2[0][1],0):
                        line2_slope = (line2[1][1]-line2[0][1])/(line2[1][0]-line2[0][0])
                        line2_perpendicular = -1/line2_slope
                        x_move = np.sqrt(1/(1+line2_perpendicular**2))
                        y_move = x_move*line2_perpendicular
                
                    elif isclose(line2[1][0]-line2[0][0],0):
                        x_move = 1
                        y_move = 0

                    else:
                        x_move = 0
                        y_move = 1

                    line2_temp = np.array(line2)

                    # Apply first move:
                    line2_temp[:,0] += x_move*nozzle_radius
                    line2_temp[:,1] += y_move*nozzle_radius

                    check_3 = self.lineCheck(line1, line2_temp, recursive_stop = True)
                    #if not self.lineCheck(line1, line2_temp, recursive_stop = True):
                    #return False

                    # Apply second move:
                    line2_temp[:,0] -= 2*x_move*nozzle_radius
                    line2_temp[:,1] -= 2*y_move*nozzle_radius

                    check_4 = self.lineCheck(line1, line2_temp, recursive_stop = True)
                    #if not self.lineCheck(line1, line2_temp, recursive_stop = True):
                    #return False

                    ###
                    # Getting an issue with a 'circular' error where one line blocks another and vice-versa

                    if not (check_1 and check_2 and check_3 and check_4) and recursive_stop_2:
                        opposite_order = self.lineCheck(line2, line1, recursive_stop = False, recursive_stop_2 = False)
                        if not opposite_order:
                            #print("Circular error")
                            return True
                        else:
                            return False

                    elif not (check_1 and check_2 and check_3 and check_4):
                        return False
                    
                    else:
                        return True
                
                else:
                    return True

    def genEdges(self,debug):
        bad_angle_lines = []
        if debug:
            logfile = open('log.csv','w')
        def edgeCheck(index,edge):
            edges_to_compare = copy.copy(self.E)
            edges_to_compare.pop(index)
            for other_edge in edges_to_compare:
                    # Check if lines intersect in the XY plane
                    # Note there was a weird numpy erray error, fixed using .tolist()
                    line_check_result = self.lineCheck([self.vert_vals[edge[0]].tolist(),
                    self.vert_vals[edge[1]].tolist()],
                    [self.vert_vals[other_edge[0]].tolist(),
                    self.vert_vals[other_edge[1]].tolist()]) 
                    if not line_check_result:
                        return (False,other_edge)
            return (True,None)

        pool = Pool(nodes=int(self.processes))
        edges = np.copy(self.E)
        index = np.arange(len(edges))
        results = pool.map(edgeCheck, index, edges)
        
        #current_edge = 1
        for index,((result,other_edge),edge) in enumerate(zip(results,edges)):
            if result:
                self.kill_list.append(index)
                if self.isValidAngle(edge[0],edge[1]) and self.isValidAngle(edge[1],edge[0]):
                    self.addEdge(edge[0],edge[1])
                else:
                    if self.isValidAngle(edge[0],edge[1]): 
                        bad_angle_lines.append([edge[0],edge[1]])
                    else:
                        bad_angle_lines.append([edge[1],edge[0]])
            else:
                if debug:
                    logfile.write('{},{},{},{},{}\n'.format(edge,other_edge,False,[self.vert_vals[edge[0]].tolist(),self.vert_vals[edge[1]].tolist()],[self.vert_vals[other_edge[0]].tolist(),self.vert_vals[other_edge[1]].tolist()]))
            #print("{}/{}".format(current_edge,len(self.E)))
            #current_edge += 1

        #pool.terminate()
        if debug:
            logfile.close()

        return bad_angle_lines

def run(edges,vertices,processes,export,filename):
    # Load mesh data from Rhino
    #F = np.loadtxt('sphere_filled_edges.csv',delimiter=',',dtype='int')
    #V = np.loadtxt('sphere_filled_vertices.csv',delimiter=',')
    F = edges
    V = vertices
    print_path = []
    print_path_temp = []

    # Create graph and fill it all the triangles from mesh
    g3 = Graph (len(V),V)
    g3.processes = processes

    # Load edges from rhino
    for f in F:
        g3.E.append(f)

    ####
    # A temp solution to the line angle problem could be to print those valid edges prior to calculating euler path
    # To do this, although they are valid edges, we would remove them from the edge generation function, printing those first or last

    debug_mode = False
    while g3.E:
    #for i in range(20):
        for kill in sorted(g3.kill_list,reverse=True):
            g3.E.pop(kill)
        g3.kill_list = []
        print('Remaining edges: {}'.format(len(g3.E)))
        bad_angle_lines = g3.genEdges(debug= debug_mode)

        # Count number of odd nodes to determine number of parts
        odd_node_count = 0
        for index in range(len(V)):
            if len(g3.graph[index])%2:
                odd_node_count += 1

        # Using a color map, create a color for each part
        colors = matplotlib.cm.rainbow(np.linspace(0, 1, odd_node_count/2))

        # Step through each color (part) and connect
        for color in colors:
            odd_nodes = []
            for index in range(len(V)):
                if len(g3.graph[index])%2:
                    odd_nodes.append(index)
            g3.color = color
            ### Issue where continuously adding odd node while failing due to print angle
            # Quick fix is to simply skip that node
            # Better fix is to update the graph
            test = g3.printEulerTour(odd_nodes[0])
            if len(test) == 1:
                import ipdb; ipdb.set_trace()
            print_path_temp.append(test)

        for index in range(len(V)):
            if len(g3.graph[index]) != 0:
                g3.color = 'k'
                # Check if print path is empty []
                if print_path_temp == []:
                    # Never seems to go through here?
                    test = g3.printEulerTour(index)
                    if len(test) == 1:
                        import ipdb; ipdb.set_trace()
                    print_path_temp.append(test)
                    #print_path_temp.append(g3.printEulerTour(index))
                else:
                    # Cycle through each of the parts
                    for part in print_path_temp:
                        # If this index is in the part
                        if index in part:
                            pos_in_part = part.index(index)
                            # Place loop in place of index
                            part[pos_in_part:pos_in_part+1] = g3.printEulerTour(index)
                            break
                        else:
                        #Cannot be appended to existing path
                            #import ipdb; ipdb.set_trace()
                            test = g3.printEulerTour(index)
                            if len(test) == 1:
                                import ipdb; ipdb.set_trace()
                            print_path_temp.append(test)
                            #print_path_temp.append(g3.printEulerTour(index))
                            break

        if not print_path_temp:
            if not debug_mode:
                debug_mode = True
            else:
                break

        #ADD BAD ANGLE LINES
        print_path_temp += bad_angle_lines
        print_path.append(print_path_temp)
        print_path_temp = []

    #### RESULTS ####

    flat_print_path = [item for sublist in print_path for item in sublist]
    flatter_print_path = [item for sublist in flat_print_path for item in sublist]

    if plot_results:
        import matplotlib.pyplot as plt
        list_length = [len(part) for part in flat_print_path]
        dist = OrderedDict(sorted(Counter(list_length).items()))
        plt.bar(range(len(dist)), list(dist.values()), align='center')
        plt.xticks(range(len(dist)), list(dist.keys()))
        plt.xlabel('Number of vertices')
        plt.ylabel('Frequency')
        plt.title('Number of vertices along continuous printed path')
        plt.savefig('output/{}_edge_distribution.png'.format(filename))

    print("Input edge count: {} Solution edge count: {}".format(len(F),len(flatter_print_path)-len(flat_print_path)))
    print("Bad angle count: {}".format(g3.bad_angle_count))

    if export:
        np.save('output/{}_path.npy'.format(filename),flat_print_path)    

    #import ipdb; ipdb.set_trace()

    return flat_print_path

if __name__ == '__main__':
    edges = np.loadtxt('sphere_filled_edges.csv',delimiter=',',dtype='int')
    vertices = np.loadtxt('sphere_filled_vertices.csv',delimiter=',')
    run(edges,vertices,sys.argv[1])