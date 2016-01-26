"""
prims_algorithm_adjacency_matrix.py

Author: Michael Munson
Date: January 25, 2016

Input: A file containing white-space separated real numbers, which are interpreted
as entries in an adjacency matrix containing distances between nodes.

Uses Prim's Algorithm to find the Minimum Spanning Forest, implemented
with a Fibonacci-heap-backed priority-queue.

Select a start node arbitrarily, and add it to the first MS Tree. At each step,
add to the MST the node nearest to another node already in the MST.

Returns a list of MST's, one for each connected component in the input
graph.


From https://en.wikipedia.org/wiki/Prim%27s_algorithm:
"For graphs that are sufficiently dense, Prim's algorithm
can be made to run in linear time, meeting or improving
time bounds for other algorithms."

Timing:

For Minimum Spanning Tree in a complete graph of 128 vertices:
~0.0067 seconds
"""

import sys
sys.path.append('./fibonacci-heap-mod-0.99')
from dijkstras_algorithm import read_cities
import numpy as np
import random
import fibonacci_heap_mod as fib
import time


# TODO: implement with adjacency list
def prims_algorithm_adjacency_matrix(adj_matrix, num_vertices):
    """
    adj_matrix is a numpy ndarray
    return minimum spanning tree (list of edges, with weights)
    """
    matrix_size = np.shape(adj_matrix)[0]

    root_index = random.randint(0, matrix_size - 1)
    mst = MST(root_index)
    num_connected_components = 1

    inf_priority = 10 ** 10 - 1

    not_included_set = set(range(matrix_size)) - {root_index}

    # set root at distance 0, all others at distance 'inf_priority'
    # heap entry is (node_index, nearest_index); entry priority is distance to nearest_index node
    not_included_heap, node_dict = initialize_not_included_heap(adj_matrix, root_index)

    # associate distance to root_index to all
    update_distances(adj_matrix, not_included_set, not_included_heap, node_dict, root_index)


    while not_included_set:
        nearest_entry = not_included_heap.dequeue_min()

        # occurs when there are multiple connected components
        if nearest_entry.m_priority == inf_priority:
            num_connected_components += 1
            root_index = random.choice(not_included_set)
            mst.roots.append(root_index)
        else:
            mst.total_dist += nearest_entry.m_priority

        new_node_index, nearest_node_index = nearest_entry.m_elem
        mst.edges[num_connected_components - 1].append((new_node_index, nearest_node_index, nearest_entry.m_priority))
        not_included_set -= {new_node_index}
        update_distances(adj_matrix, not_included_set, not_included_heap, node_dict, new_node_index)

    return mst

def initialize_not_included_heap(adj_matrix, root_index):
    not_included_heap = fib.FibonacciHeap()
    node_dict = {}

    inf_priority = 10 ** 10 - 1

    for i in range(0, np.shape(adj_matrix)[0]):
        # element is: (distance to nearest node in current tree, index of this node)
        if i != root_index:
            node_dict[i] = not_included_heap.enqueue((i, root_index), inf_priority)

    return not_included_heap, node_dict


# TODO: deal with possibility of inf distance
def update_distances(adj_matrix, not_included_set, not_included_heap, node_dict, node_index):
    distances = adj_matrix[node_index, :]

    for i in not_included_set:
        new_distance = distances[i]
        old_min_distance = node_dict[i].m_priority

        if new_distance < old_min_distance:
            node_dict[i].set_value((i, node_index))  # node at node_index is now nearest node to i
            not_included_heap.decrease_key(node_dict[i], new_distance)
    return


class MST(object):
    def __init__(self, root_index):
        self.roots = [root_index]
        self.edges = [[]]
        self.total_dist = 0

cities = read_cities()

src = 0
num_cities = 128

# #Timing code:
# min_elapsed = 10 ** 10 - 1
# for _ in range(10):
#     start = time.clock()
#     results_list = prims_algorithm_adjacency_matrix(cities, num_cities)
#     end = time.clock()
#     elapsed = end - start
#     if elapsed < min_elapsed: min_elapsed = elapsed
#
# print ("Min elapsed time is {0}".format(min_elapsed))

mst = prims_algorithm_adjacency_matrix(cities, num_cities)

print("Total of {0} MSTs, rooted at {1}\nTotal distance is {2} miles\nEdges are:\n{3}".
      format(len(mst.roots), mst.roots,mst.total_dist, mst.edges))
