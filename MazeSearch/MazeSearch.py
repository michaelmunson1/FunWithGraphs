
def a_star_search(m_graph):
    path = []

    return path


def maze_search(width, height, walls, start, end):

    walls_dict = walls_to_dict(walls)
    m_graph = MazeGraph(width, height, walls_dict, start, end)

    # path = a_star_search(m_graph)

    #return path

    # Print to check correct initialization of graphs:
    #
    # for i in range(width):
    #     for j in range(height):
    #         if (i, j) not in walls:
    #             print "adjacencies to (" + str(i) + "," + str(j) + "):"
    #             print m_graph.adjs[(i, j)]

    return m_graph


class MazeGraph:

    def __init__(self, width, height, walls_dict, start, end):  # x,y coordinates of walls, start, & end are 0-indexed

        self.width = width
        self.height = height
        self.start = start
        self.end = end
        self.adjs = {}      # TODO:

        for i in range(height):       # row
            for j in range(width):    # col
                if (i, j) in walls_dict:
                    print("got here")
                    continue
                else:
                    this_adj = []

                    deltas = [[0, 1], [1, 0], [0, -1], [-1, 0]]    # (x,y) deltas for surrounding squares (N,E,S,W)

                    for d in deltas:
                        new_col = j+d[0]
                        new_row = i+d[1]
                        if new_col >= 0 and new_col < width and new_row >= 0 and new_row < height \
                                and (new_col, new_row) not in walls_dict:

                            this_adj.append((new_col, new_row))
                    self.adjs[(j, i)] = this_adj


class SearchNode:

    def __init__(self):
        return


def walls_to_dict(walls):
    w_dict = {}

    for w in walls:
        w_dict[w] = 1

    return w_dict


def main():
    # TODO: Make test loop here, generate random input values and test

    width = 3
    height = 3
    walls = [(1, 1)]
    start = [(0, 0)]
    end = [(2, 2)]

    path = maze_search(width, height, walls, start, end)

    # print path


if __name__ == "__main__":
    main()