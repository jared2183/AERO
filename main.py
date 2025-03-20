# Python program print Eulerian Trail in a given Eulerian or Semi-Eulerian Graph
# Required packages:
# * Pathos
# * Numpy
# * Matplotlib
# * Vpython

from collections import defaultdict, OrderedDict, Counter
import copy
import numpy as np
from pathos.multiprocessing import ProcessingPool as Pool
import sys
from mecode.main import G
import pathPlanner

def gcodeGeneration(vertices, print_path, filename, travel_height=30,
                    travel_feed=2, print_feed=1, com_port=6, title="AERO",
                    volumetric_mode=False, flowrate=0.007, preview=False):
    """
        Parameters
        ----------
        vertices : list of list of floats
            A list of coordinates for the vertices [(x_1,y_1,z_1),(x_2,y_2,z_2),...]
        print_path : list of list of ints
            A list of the order of the print path [(v_1,v_2,v_3,v_4),(v_5,v_6,...),...]
        filename : path or None (default: None)
            If a path is specified, the compiled gcode will be written to that
            file.
        travel_height : float
            Safe plane in which to travel between printing moves.
        travel_feed : float
            Speed at which to move when travelling.
        print_feed : float
            Speed at which to move when printing.
        com_port : int or None (default: 6)
            Serial communication port to use when toggle pressure supply.
        title : str or None (default: "AERO")
            Optional title added as comment to start of gcode
        volumetric_mode : bool (default: False)
            When True, will generate an extra axis command to control a Viscotec volumetric extruder
        volumetric_axis : str (default: "a")
            Which axis to use for the volumetric extruder
        flowrate : float (default: 0.007)
            Flowrate in ml/s to control the volumetric extruder  
    """
    
    #Generate gcode for path
    g = G(outfile='{}.pgm'.format(filename),print_lines=False,aerotech_include=True)
    V = vertices

    #Path settings
    end_extra_length = 0.0
    repeated_extra_height = 0.0

    #Path variables
    volumetric_extrusion = 0
    visited_vertices = []

    g.write("; "+ title)
    if volumetric_mode:
        #Zero extruder if not already done.
        g.write('G92 a0.0')
    g.write("VELOCITY OFF")
    g.rename_axis(z='A')
    g.absolute()
    g.feed(travel_feed)

    def getPosition():
        return np.fromiter(g.current_position.values(), dtype=float)[:3]

    def nodeDelay(volumetric_extrusion,delay_time=1.5):
        volumetric_extrusion += flowrate*delay_time
        if volumetric_mode:
            g.write('G1 a{:.6f} E{:.6f}'.format(volumetric_extrusion,flowrate))
        else:
            g.dwell(delay_time)
        return volumetric_extrusion

    def startDelay(volumetric_extrusion,delay_time=1.5):
        volumetric_extrusion += flowrate*delay_time
        if volumetric_mode:
            g.write('G1 a{:.6f} E{:.6f}'.format(volumetric_extrusion,flowrate))
        else:
            g.dwell(delay_time)
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

    #### NOTES
    # Node intersections seem to cause seperation in X and Y, but not so much
    # in z. This would be caused by the translating nozzle into the node and
    # less by the nozzle pulling away from the node. Solution: remove z offset, 
    # just keep x and y
    # Pretty much elongate end of along line
    # and push out nodes in combination of line extension
    # Seems to be working correctly thus far

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

                # checks if there should be an offset
                if print_line.any() and end_extra_length != 0:
                    mag_line = np.sqrt(print_line[0]**2+print_line[1]**2+print_line[2]**2)
                    unit_vector = print_line/mag_line
                    offset_vertex = np.array(V[edge])+unit_vector*end_extra_length
                else:
                    offset_vertex = np.array(V[edge])
                
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
    if preview:
        #g.view('matplotlib',color_on=False)
        g.view('vpython',nozzle_cam=False)
    g.teardown()

if __name__ == '__main__':
    #import files
    model = 'octet'
    edges = np.loadtxt('examples/Octet Lattice 3x3x3/Input/octet_edges.csv',delimiter=',',dtype='int')
    vertices = np.loadtxt('examples/Octet Lattice 3x3x3/Input/octet_vertices.csv',delimiter=',')
    # Run Path Planning
    parallel_nodes = 8
    path = pathPlanner.run(edges,vertices,processes=parallel_nodes,export=True,filename=model)

    # Or input existing path
    #path = np.load('octet_path.npy',allow_pickle=True)

    gcode = gcodeGeneration(vertices,path,filename=model)