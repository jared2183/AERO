import numpy as np
from main_export import G

#Generate gcode for path
travel_height = 53.0
travel_feed = 20.0
print_feed = 0.5
com_port = 6
g = G(outfile='sphere_lattice_filled_v3.pgm',print_lines=False,aerotech_include=True)
F = np.loadtxt('Bunny/sphere_filled_edges.csv',delimiter=',',dtype='int')
V = np.loadtxt('Bunny/sphere_filled_vertices.csv',delimiter=',')
print_path = np.load('sphere_filled_print_path_Final.npy')

#F = np.load('edges_FCC.npy')
#V = np.load('vertex_FCC.npy')
#print_path = np.load('print_path_FCC.npy') 

#end_extra_length = 0.8
#node_extra_length = 0.6
end_extra_length = 0.0
node_extra_length = 0.0
repeated_extra_height = 0.0
volumetric_extrusion = 0
visited_vertices = []
flowrate = 0.004
volumetric_mode = False

g.write("; Bunny Lattice V1.36 by rweeks on 01/29/19")
g.write("VELOCITY OFF")
#Zero extruder if not already done
if volumetric_mode:
    g.write('G92 a0.0')
g.rename_axis(z='A')
g.absolute()
g.feed(travel_feed)

def getPosition():
    return np.fromiter(g.current_position.values(), dtype=float)[:3]

def nodeDelay(volumetric_extrusion,delay_time=1.5):
    # 2 second dwell
    volumetric_extrusion += flowrate*delay_time
    if volumetric_mode:
        g.write('G1 a{:.6f} E{:.6f}'.format(volumetric_extrusion,flowrate))
    else:
        g.dwell(2)
    return volumetric_extrusion

def startDelay(volumetric_extrusion,delay_time=2.0):
    # 4 second dwell
    volumetric_extrusion += flowrate*delay_time
    if volumetric_mode:
        g.write('G1 a{:.6f} E{:.6f}'.format(volumetric_extrusion,flowrate))
    else:
        g.dwell(4)
    return volumetric_extrusion

def suckBack(volumetric_extrusion):
    volumetric_extrusion -= 0.05
    if volumetric_mode:
        g.write('G1 a{:.6f} E0.05'.format(volumetric_extrusion))
    return volumetric_extrusion

def antiSuckBack(volumetric_extrusion):
    volumetric_extrusion += 0.03
    if volumetric_mode:
        g.write('G1 a{:.6f} E0.05'.format(volumetric_extrusion))
    return volumetric_extrusion


### Want to apply offsets
# Pretty much elongate end of along line
# and push out nodes in combination of line extension
# Seems to be working correctly thus far

#### NOTES
# Node intersections seem to cause seperation in X and Y, but not so much
# in z. This would be caused by the translating nozzle into the node and
# less by the nozzle pulling away from the node. Solution: remove z offset, 
# just keep x and y

# Decrease flowrate a bit, starting to look kind of messy compared to original bunny

# End points need to still be extended a bit further, could likely just keep it in x and y as well

# Changes:
# Decreased flowrate by 20%
# Decreased start delay by 25%
# Increased end_extra_length by %20
# Nodes offsets only occur in X&Y and not Z
# End offsets only occur in X&Y and not Z


for part in print_path:
    g.write("G65 F100")
    g.move(*V[part[0]][:2])
    g.write("G66 F1")
    if part[0] in visited_vertices:
        V[part[0]][2] += repeated_extra_height
    else:
        visited_vertices.append(part[0]) 
    
    g.move(z=V[part[0]][2])
    g.write("G65 F1")
    g.feed(print_feed)
    g.toggle_pressure(com_port)
    volumetric_extrusion =startDelay(volumetric_extrusion)
    for index,edge in enumerate(part[1:]):
        if index == len(part[1:]) - 1:
            #On last move
            #For the last move, we simply want to extend along line
            # Get unit vector for line and add to final position
            print_line = np.array(V[edge]-getPosition())
            # #Remove vertical component
            print_line[2] = 0
            mag_line = np.sqrt(print_line[0]**2+print_line[1]**2+print_line[2]**2)
            unit_vector = print_line/mag_line
            offset_vertex = np.array(V[edge])+unit_vector*end_extra_length
            
            if edge in visited_vertices:
                offset_vertex[2] += repeated_extra_height
            else:
                visited_vertices.append(edge)

            offset_vertex_line = offset_vertex-np.array(getPosition())
            time_for_move = np.sqrt(offset_vertex_line[0]**2+offset_vertex_line[1]**2+offset_vertex_line[2]**2)/print_feed
            volumetric_extrusion += time_for_move*0.004
            if volumetric_mode:
                g.move(*offset_vertex,a=volumetric_extrusion)
            else:
                g.move(*offset_vertex)
            volumetric_extrusion = nodeDelay(volumetric_extrusion)
            volumetric_extrusion = suckBack(volumetric_extrusion)

        else:
            #Going into node
            # Get unit vector for line entering node
            print_line_in = np.array(V[edge]-getPosition())
            # Remove vertical component
            print_line_in[2] = 0
            mag_line_in = np.sqrt(print_line_in[0]**2+print_line_in[1]**2+print_line_in[2]**2)
            mag_line_in = 1
            unit_vector_in = print_line_in/mag_line_in

            # Get unit vector for line exiting node
            print_line_out = np.array(V[part[1:][index+1]]-V[edge])
            # Remove vertical component
            print_line_out[2] = 0
            mag_line_out = np.sqrt(print_line_out[0]**2+print_line_out[1]**2+print_line_out[2]**2)
            mag_line_out = 1
            unit_vector_out = print_line_out/mag_line_out

            offset_unit_vector = unit_vector_in - unit_vector_out
            offset_vertex = np.array(V[edge])+offset_unit_vector*end_extra_length

            if edge in visited_vertices:
                offset_vertex[2] += repeated_extra_height
            else:
                visited_vertices.append(edge)

            offset_vertex_line = offset_vertex-np.array(getPosition())
            time_for_move = np.sqrt(offset_vertex_line[0]**2+offset_vertex_line[1]**2+offset_vertex_line[2]**2)/print_feed
            volumetric_extrusion += time_for_move*0.004
            if volumetric_mode:
                g.move(*offset_vertex,a=volumetric_extrusion)
            else:
                g.move(*offset_vertex)
            volumetric_extrusion = nodeDelay(volumetric_extrusion)

    g.toggle_pressure(com_port)
    g.feed(travel_feed)
    g.write("G66 F100")
    g.move(z=travel_height)
    volumetric_extrusion =antiSuckBack(volumetric_extrusion)

g.home()
g.write("VELOCITY ON")
#g.view('matplotlib',color_on=True)
#g.view('vpython',nozzle_cam=False)
g.teardown()