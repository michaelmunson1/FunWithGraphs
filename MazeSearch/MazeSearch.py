import math
import Queue


def a_star_search(m_graph):
    """
    :rtype: path (list)
    """

    start_node = SearchNode(m_graph.start, None, 0, manhattan_distance(m_graph.start, m_graph.end))  # no parent
    path = []
    frontier = Queue.PriorityQueue()
    #frontier.put(start_node.path_cost, (start_node, start_node.path_cost))  # path_cost is 0

    frontier_dict = {}
    frontier_dict[start_node.square] = 1
    explored_set = {}

    print ("Exploring node at location " + str(m_graph.start))

    path = a_star_inner(m_graph, path, start_node, explored_set, frontier, frontier_dict)

    return path


def a_star_inner(m_graph, path, curr_node, explored_set, frontier, frontier_dict):
    """
    :rtype: path (list)
    """
    if curr_node.square == m_graph.end:

        print ("\nSolution found!: \n")
        path = print_path_to(curr_node)

        return path

    else:
        # find child node of current nearest to end (manh dist)

        explored_set[curr_node.square] = 1

        best_square = None
        best_dist = 999999999

        for next_square in m_graph.adjs[curr_node.square]:  # find square nearest to goal location
            if next_square not in frontier_dict:
                if next_square not in explored_set:
                    dist = manhattan_distance(next_square, m_graph.end)
                    path_cost = curr_node.path_cost + 1
                    total_cost = path_cost + dist  # f = g + h  (cost to this node + heuristic estimate)

                    next_node = SearchNode(next_square, curr_node, path_cost, total_cost)

                    frontier.put((total_cost, next_node))
                    frontier_dict[next_square] = 1

        if not frontier.empty():
            (total_cost, next_node) = frontier.get()
            next_square = next_node.square
            del frontier_dict[next_square]
            print ("Exploring node at location " + str(next_square))

            return a_star_inner(m_graph, path, next_node, explored_set, frontier, frontier_dict)
        else:

            return None   # no solution found

        # for next_square in m_graph.adjs[curr_node.square]:
        #     if next_square not in explored_set:
        #
        #         if dist < best_dist:
        #             best_square = next_square
        #             best_dist = dist
        #
        # if best_square is not None:  # at least one adj node not yet explored
        #     next_node = SearchNode(best_square, curr_node, curr_node.path_cost + 1)
        #     path.append(next_node.square)
        #
        #     print ("Exploring node at location " + str(best_square))
        #
        #     return a_star_inner(m_graph, path, next_node, explored_set)
        # else:  # all adj nodes explored - backtrack
        #     next_node = curr_node.parent
        #     if next_node is None:  # curr_node is the start node
        #         return None  # no solution found
        #     else:
        #         del path[-1]
        #         return a_star_inner(m_graph, path, next_node, explored_set)

def print_path_to(node):
    path_reverse = []
    path = []

    while node is not None:
        path_reverse.append(node.square)
        node = node.parent

    path_len = len(path_reverse)

    print("Length is " + str(path_len - 1) + "\n")

    for i in range(path_len):
        square = path_reverse[path_len - i - 1]
        path.append(square)
        print square

    return path


def manhattan_distance(loc1, loc2):
    return int(math.fabs(loc1[0] - loc2[0]) + math.fabs(loc1[1] - loc2[1]))


def maze_search(width, height, walls, start, end):
    """ CHECK FOR BOUNDARY CASES:
            width or height is <= 0
            start or end has negative coordinate or is not given
    """

    walls_dict = walls_to_dict(walls)
    m_graph = MazeGraph(width, height, walls_dict, start, end)

    path = a_star_search(m_graph)

    return path

    # Print to check correct initialization of graphs:
    #
    # for i in range(width):
    #     for j in range(height):
    #         if (i, j) not in walls:
    #             print "adjacencies to (" + str(i) + "," + str(j) + "):"
    #             print m_graph.adjs[(i, j)]


class MazeGraph:
    def __init__(self, width, height, walls_dict, start, end):  # x,y coordinates of walls, start, & end are 0-indexed

        self.width = width
        self.height = height
        self.start = start
        self.end = end
        self.adjs = {}  # TODO:

        for i in range(height):  # row
            for j in range(width):  # col
                if (j, i) in walls_dict:
                    continue

                this_adj = []

                deltas = [[0, 1], [1, 0], [0, -1], [-1, 0]]  # (x,y) deltas for surrounding squares (N,E,S,W)

                for d in deltas:
                    new_col = j + d[0]
                    new_row = i + d[1]
                    if new_col >= 0 and new_col < width and new_row >= 0 and new_row < height \
                            and (new_col, new_row) not in walls_dict:
                        this_adj.append((new_col, new_row))
                self.adjs[(j, i)] = this_adj


class SearchNode:
    def __init__(self, square, parent, path_cost, total_cost):
        self.square = square  # type: (x,y)
        self.parent = parent  # type: SearchNode
        self.path_cost = path_cost
        self.total_cost = total_cost


def walls_to_dict(walls):
    w_dict = {}

    for w in walls:
        w_dict[w] = 1

    return w_dict


def main():
    # TODO: Make test loop here, generate random input values and test

    width = 10
    height = 10
    walls = [(7, 0), (8, 0), (0, 1), (2, 1), (3, 1), (7, 1), (8, 1), (9, 1), (0, 2), (2, 2),
             (4, 2), (6, 2), (2, 3), (4, 3), (1, 4), (2, 4), (5, 4), (7, 4), (8, 4), (1, 5), (2, 5),
             (3, 6), (4, 6), (5, 6), (8, 6), (0, 7), (2, 7), (5, 7), (6, 7), (8, 7), (9, 7),
             (2, 8), (0, 9), (6, 9), (7, 9)]
    start = (0, 0)
    end = (9, 9)

    path = maze_search(width, height, walls, start, end)

    # print path


if __name__ == "__main__":
    main()
