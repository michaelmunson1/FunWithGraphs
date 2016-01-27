# -*- coding: utf-8 -*-
"""
kruskals_algorithm_adjacency_matrix.py

Author: Michael Munson
Date: January 25, 2016

Input: A file containing white-space separated real numbers, which are interpreted
as entries in an adjacency matrix containing distances between nodes.

Uses Kruskal's Algorithm to find the Minimum Spanning Forest of a graph.
Uses counting sort to sort edges by weight, and a disjoint-set data-structure.
The time complexity is O(E α(V)), where α is the extremely slowly growing inverse
of the single-argument Ackermann function (see
https://en.wikipedia.org/wiki/Ackermann_function#Definition_and_properties).

Begin by separating each node into its own tree. Select the lowest-weight edge
that joins nodes in distinct trees, and join those trees. Repeat until all
components are connected.

Returns a list of MST's, one for each connected component in the input
graph.

Implementation of Kruskal's algorithm is somewhat simpler than Prim's
(which requires storing and updating several storage tables), with
 essentially the same asymptotics

from https://en.wikipedia.org/wiki/Prim%27s_algorithm:
"For graphs that are sufficiently dense, Prim's algorithm
can be made to run in linear time, meeting or improving
time bounds for [Boruvka and Kruskal's] algorithms."



Timing:

For Minimum Spanning Tree in a complete graph of 128 vertices:

"""

def kruskals_algorithm(adj_matrix):
    sorted_edge_list = adj_matrix_to_edge_list(adj_matrix)

    return MSF

def adj_matrix_to_edge_list(adj_matrix):
    indices = range(np.shape(adj_matrix)[0])
    edge_list = [(i, j, adj_matrix[i, j]) for i in indices for j in indices]

    sorted_edge_list = counting_sort(edge_list)

    return sorted_edge_list

class FindUnion(object):
    """
    For storing disjoint sets
    """
    def __init__(self):
        return

    # TODO: Proper to place comparison here, or in Node?
    def __cmp__(self, other):



class Node(object):
    """
    To be used in FindUnion disjoint sets data structure
    """
    def __init__(self):
        self.parent = self
