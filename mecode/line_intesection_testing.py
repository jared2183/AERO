### Test Cases ####
#### PARALLEL in XY ####
# Parallel no overlap
"""
A = [0.0,0.0,0.0]
B = [1.0,1.0,0.0]
C = [1.5,1.5,0.0]
D = [2.0,2.0,0.0]
"""

# Parallel overlap, line2 accessible
"""
A = [0.0,0.0,0.0]
B = [1.0,1.0,0.0]
C = [0.5,0.5,0.5]
D = [2.0,2.0,0.5]
"""

# Parallel overlap, line2 inaccessible
"""
A = [0.0,0.0,0.5]
B = [1.0,1.0,0.5]
C = [0.5,0.5,0.0]
D = [2.0,2.0,0.0]
"""

# Lines are parallel in XY and share the same node with line2 accessible
"""
A = [0.0,0.0,0.0]
B = [3.5,0.0,0.0]
C = [0.0,0.0,0.0]
D = [3.5,0.0,3.5]
"""

# Lines are parallel in XY and share the same node with line2 inaccessible
"""
A = [0.0,0.0,0.0]
B = [3.5,0.0,3.0]
C = [0.0,0.0,0.0]
D = [3.5,0.0,0.0]
"""

# Lines are parallel in XY and share the same node but face opposing directions
"""
A = [0.0,0.0,0.0]
B = [3.5,0.0,3.5]
C = [3.5,0.0,3.5]
D = [7.0,0.0,0.0]
"""

# Lines are parallel and overlap, but one of the lines is vertical

"""
A = [0.0,0.0,0.0]
B = [0.0,0.0,3.5]
C = [0.0,0.0,0.0]
D = [0.0,3.5,3.5]
"""

#### NOT PARALLEL in XY ####
# Line intersect and segment intersect, line2 passes over line1
"""
A = [-0.5,-0.5,0.3]
B = [1.0,1.0,0.0]
C = [1.0,1.0,1.0]
D = [-1.5,3.5,1.3]
"""

# Line intersect and segment intersect, line2 passes under line1
"""
A = [-0.5,-0.5,1.3]
B = [1.0,1.0,1.0]
C = [1.0,1.0,0.0]
D = [-1.5,3.5,0.3]
"""

# Line intersect and segment intersect and collide3D
"""
A = [-0.5,-0.5,0.3]
B = [1.0,1.0,1.0]
C = [1.0,1.0,1.0]
D = [-1.5,3.5,1.3]
"""

# Line intersect not segment intersect
"""
A = [-0.5,-0.5,0.0]
B = [1.0,1.0,0.0]
C = [1.0,2.5,0.0]
D = [-1.5,3.5,0.0]
"""

# Line intersect and segment intersect with shared node in 3D
"""
A = [0.0,0.0,0.0]
B = [3.5,0.0,0.0]
C = [0.0,0.0,0.0]
D = [0.0,3.5,3.5]
"""

# Line intersect and segment intersect with shared node in 3D
"""
A = [0.0,0.0,0.0]
B = [3.5,0.0,3.5]
C = [0.0,0.0,7.0]
D = [0.0,3.5,3.5]
"""

# Line intersect, but segment doesn't intersect second line is vertical
"""
A = [0.0,0.0,0.0]
B = [3.5,0.0,3.5]
C = [0.0,7.0,0.0]
D = [0.0,3.5,3.5]
"""

# Other

A = [7.0, -7.0, -7.0]
B = [7.0, -7.0, 0.0]
C = [7.0, -7.0, -7.0]
D = [7.0, 0.0, -7.0]



import matplotlib.pyplot as plt
import numpy as np

plt.plot([A[0],B[0]],[A[1],B[1]])
plt.plot([C[0],D[0]],[C[1],D[1]])
plt.show()

def lineCheck(line1, line2):
    #Check if printing line1 will prevent printing line2
    def onSegment(p,q,r):
        #Check if point r is on line segment pq
        def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
        #Function for comparing floats
        #https://stackoverflow.com/questions/5595425/what-is-the-best-way-to-compare-floats-for-almost-equality-in-python
            return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
        
        try:
            line_slope = (q[1]-p[1])/(q[0]-p[0])
            constant = p[1]-line_slope*p[0]
            
            if isclose(r[1]-line_slope*r[0],constant):
                #Point is on line
                #Check if point is in segment
                if (r[0] <= max(p[0], q[0]) and r[0] >= min(p[0], q[0]) and
                    r[1] <= max(p[1], q[1]) and r[1] >= min(p[1], q[1])):
                    return True
        except ZeroDivisionError:
            #Line being compared is vertical (infinite slope)
            #Check if point is in segment
            if r[0] == p[0] and r[1] <= max(p[1], q[1]) and r[1] >= min(p[1], q[1]):
                return True
        return False

    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)

    ### PARALLEL CASE ###
    if div == 0: #Lines are parallel
        #Check if line1 overlaps line2, including endpoints
        # if 0 endpoints are on segment they must not overlap
        # If one endpoint is on segment they must be opposing
        # if two endpoints are on segment they must be overlapping (or could be a vertical line with another line)
        if (onSegment(line2[0],line2[1],line1[0]) and onSegment(line2[0],line2[1],line1[1])) or (onSegment(line1[0],line1[1],line2[0]) and onSegment(line1[0],line1[1],line2[1])):
            # Line1 overlaps line2
            # Check if one of the lines is vertical
            def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
            #Function for comparing floats
            #https://stackoverflow.com/questions/5595425/what-is-the-best-way-to-compare-floats-for-almost-equality-in-python
                return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
            #import ipdb; ipdb.set_trace()
            if (isclose(line1[0][0],line1[1][0]) and isclose(line1[0][1],line1[1][1])) or (isclose(line2[0][0],line2[1][0]) and isclose(line2[0][1],line2[1][1])):
                #One of the lines is vertical
                print 'lines are parallel and overlap, but one of the lines is vertical'
                # GOOD to print
                return True
            else:
                # No vertical line
                # Solution assumes no intersection other plane when overlapping in XY plane
                # GOOD TO PRINT LINE1 IF LOWER IN Z THAN LINE2 (LINE 2 STILL ACCESSIBLE)
                # Since they are parallel as viewed from XY plane (in same plane), whichever average Z is more will be higher
                average_line1 = (line1[0][2] + line1[1][2])/2.0
                average_line2 = (line2[0][2] + line2[1][2])/2.0
                if average_line2 > average_line1:
                    print 'lines are parallel and overlap, but line2 still ACCESSIBLE'
                    # GOOD to print
                    return True
                else:
                    print 'lines are parallel and overlap and line2 is NOT ACCESSIBLE'
                    # NOT GOOD to print
                    return False
        elif onSegment(line2[0],line2[1],line1[0]) or onSegment(line2[0],line2[1],line1[1]):
            # Line1 overlaps line2
            print "lines are parallel but facing opposing directions"
            # GOOD to print
            return True
        else:
            # Lines are parallel but don't overlap
            # PRINTING LINE1 WILL NOT STOP LINE2 FROM BEING PRINTED
            print 'lines are parallel and dont overlap'
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
            unique_points = np.unique(all_points,axis=0)
            num_unique_points = len(unique_points)
            if num_unique_points < num_points: #There must have been duplicates
                #Points share a node
                print 'lines share the same node, but arent parallel, thus printing line1 cant block line2'
                # GOOD to print
                return True
            # GOOD TO PRINT LINE1 IF LOWER IN Z THAN LINE2 (LINE 2 STILL ACCESSIBLE)
            # Cannot just use average, since not in same plane
            # Line2 is higher than line1 if Z value is higher at XY coordinate of intersection
            fraction_of_line_at_intersect_line1 = abs(intersect_point[0]-line1[0][0])/abs(line1[1][0]-line1[0][0])
            fraction_of_line_at_intersect_line2 = abs(intersect_point[0]-line2[0][0])/abs(line2[1][0]-line2[0][0])
            height_at_intersect_line1 = fraction_of_line_at_intersect_line1*(line1[1][2]-line1[0][2])+line1[0][2]
            height_at_intersect_line2 = fraction_of_line_at_intersect_line2*(line2[1][2]-line2[0][2])+line2[0][2]
            def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
                #Function for comparing floats
                #https://stackoverflow.com/questions/5595425/what-is-the-best-way-to-compare-floats-for-almost-equality-in-python
                return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
            if height_at_intersect_line2 > height_at_intersect_line1:
                #Line2 passes over line1, thus it is fine to print line1
                # GOOD to print
                print 'segments intersect, but line2 passes over line1'
                return True
            elif isclose(height_at_intersect_line1,height_at_intersect_line2):
                #Lines either share the same node or could have collied in mid-air
                #Either way it should be fine to print, no benefit either way
                # GOOD to print
                print 'segments intersect, lines actually collided in 3d'
                return True
            else:
                #Line2 passes under line1, thus printing line1 would block access to line2
                # NOT GOOD to print
                print 'segments intersect, but line2 passes under line1'
                return False
        else:
            # Segments don't intersect
            print "segments dont intersect"
            # GOOD to print
            return True

#plt.plot([A[0],B[0]],[A[1],B[1]])
#plt.plot([C[0],D[0]],[C[1],D[1]])
print lineCheck((A, B), (C, D))
#plt.show()

