"""
dijkstras_algorithm.py

Author: Michael Munson
Date: January 24, 2016

Uses pairing heap

Timing:

For shortest path for 127 destinations, in a complete graph of 128 vertices:
~0.0157 seconds, vs. ~0.0089 seconds for a fibonacci-heap implementation.


"""

# TODO: timing run vs. Fibonacci heap?

import numpy as np
import pairing_heap as ph
import time


class RandomUndirectedGraph(object):
    def __init__(self, num_vertices, density):
        assert density < 1


class NodeEntry(object):
    def __init__(self, dist=float('inf'), visited=False, prev_index=None, num_steps=0):
        self.dist = dist
        self.visited = visited
        self.prev_index = prev_index
        self.num_steps = num_steps


def read_cities():
    cities_list = []
    with open('../city_distances.txt') as city_dist_file:

        for line in city_dist_file:
            line = line.strip()

            if line[0] != '#':
                this_city= [int(d) for d in line.split()]
                # print(distances)
                cities_list.append(this_city)
    cities = np.array(cities_list)
    return cities


# TODO: finish - in order to find all shortest distances
def no_finite_dist_unvisited_cities(results_list, unvisited):
    not_found_finite_dist_city = True

    for i in unvisited:
        if results_list[i].dist < float('inf'):
            not_found_finite_dist_city = False
            break

    return not_found_finite_dist_city


def dijkstras_algorithm(cities, num_cities, src_index, dest_index=None):
    """
    :param cities: numpy array of distances between cities
    :param num_cities:
    :param src_index:
    :param dest_index: When 'None', find shortest route to all cities
    :return: results_list
    """
    # TODO: describe results_list

    curr_src_index = src_index

    unvisited_node_dict = {}

    results_list = [NodeEntry() for _ in range(num_cities)]
    results_list[src_index].dist = 0
    results_list[src_index].visited = True     # TODO: needed, for any nodes?
    unvisited = set(range(num_cities)) - {src_index}

    unvisited_heap = ph.PairingHeap()

    while unvisited:

        min_dist = float('inf')
        min_dist_index = None

        #print("src index is {0} and dist is {1}".format(curr_src_index, results_list[curr_src_index].dist))

        base_distance = results_list[curr_src_index].dist

        for i in unvisited:          # TODO: in general case, use unvisited neighbors (intersect sets)

            # find minimum total distance to city i through curr_src_index
            dist_curr_to_i = cities[curr_src_index][i]
            total_distance = base_distance + dist_curr_to_i

            # update distance and route to city i
            if total_distance < results_list[i].dist:
                results_list[i].dist = total_distance

                if i not in unvisited_node_dict:
                    # store node in dict for later call to adjust_key()
                    unvisited_node_dict[i] = unvisited_heap.insert((total_distance, i))
                else:
                    unvisited_heap.adjust_key(unvisited_node_dict[i], (total_distance, i))

                results_list[i].prev_index = curr_src_index
                results_list[i].num_steps = results_list[curr_src_index].num_steps + 1

        # # find unvisited city with least distance
        # for city_index in unvisited:
        #     this_dist = results_list[city_index].dist
        #     if this_dist < min_dist:
        #         min_dist = this_dist
        #         min_dist_index = city_index

        min_dist, min_dist_index = unvisited_heap.extract()

        # update next node to 'visited'
        results_list[min_dist_index].visited = True    # TODO: needed?
        unvisited -= set([min_dist_index])
        del unvisited_node_dict[min_dist_index]

        # done
        if no_finite_dist_unvisited_cities(results_list, unvisited) or dest_index and results_list[dest_index].visited:
            return results_list

        curr_src_index = min_dist_index


def recover_path(results_list, src, dest):
    path = []
    city_index = dest
    distance = distance = results_list[dest].dist
    path.append((city_index, distance))

    while city_index:
        prev_index = results_list[city_index].prev_index
        if prev_index:
            # print(prev_index)
            path.append((prev_index, results_list[prev_index].dist))
        city_index = prev_index

    path.append((city_index, results_list[city_index].dist))

    path.reverse()
    return path


def check_path(path, cities, src, dest):
    checks_out = True

    if src != path[0][0] or path[0][1] != 0:
        #print("YES! Correct start")
    #else:
        #print("NO! Incorrect start")
        checks_out = False

    prev_city = src
    prev_city_dist = 0

    for step in path[1:]:
        curr_city = step[0]
        curr_dist = step[1]

        if prev_city_dist + cities[prev_city][curr_city] != curr_dist:
            #print("YES! Step from {0} to {1} is correct".format(prev_city, curr_city))
        #else:
            #print("NO! Step from {0} to {1} is incorrect".format(prev_city, curr_city))
            checks_out = False

        prev_city = curr_city
        prev_city_dist = curr_dist

    if path[-1][0] != dest:
        #print("YES! Correct dest")
    #else:
        #print("NO! Incorrect dest")
        checks_out = False

    return checks_out


cities = read_cities()

src = 0
num_cities = 128

# dest not set - will find shortest distance to all cities from src city

for _ in range(10):
    start = time.clock()
    results_list = dijkstras_algorithm(cities, num_cities, src)
    end = time.clock()
    print ("total time elapsed is {0}".format(end-start))

results_list = dijkstras_algorithm(cities, num_cities, src)

max_steps = max_dist = max_steps_index = max_dist_index = -1

for i in range(num_cities):
    if results_list[i].num_steps > max_steps:
        max_steps = results_list[i].num_steps
        max_steps_index = i
    if results_list[i].dist > max_dist:
        max_dist = results_list[i].dist
        max_dist_index = i

print('\nmax_steps is {0} for city # {1}'.format(max_steps, max_steps_index))
print('max_dist is {0} for city # {1}\n'.format(max_dist, max_dist_index))

path = recover_path(results_list, src, max_steps_index)
print("Printing path with maximum number of steps:\n\n{0}".format(path))
print("Path is correct: {0}\n".format(check_path(path, cities, src, max_steps_index)))

path = recover_path(results_list, src, max_dist_index)
print("Printing path with maximum total distance:\n\n{0}".format(path))
print("Path is correct: {0}".format(check_path(path, cities, src, max_dist_index)))