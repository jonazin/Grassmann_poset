import numpy as np
import itertools
import networkx as nx

def gauss_binom(n,k,q):
    numerator = np.prod([((q ** (n-j))-1) for j in range(k)])
    denominator = np.prod([((q ** (j+1))-1) for j in range(k)])
    return int(numerator / denominator)


def build_grassmann_poset(n, d, q):
    F = GF(q)
    V = VectorSpace(GF(q),n)
    P = set([])
    P_by_dims = {}
    for i in range(d+1):
        P_by_dims[i] = set(V.subspaces(i))
        P = P.union(P_by_dims[i])

#     expected_size = np.sum([gauss_binom(n,k,q) for k in range(d+1)])
#     print("the expected poset size is {}".format(expected_size))
#     for t in itertools.product(V,repeat=d):
#         U = V.subspace([t[i] for i in range(d)])
#         P_by_dims[U.dimension()].add(U)
#         if all([len(P_by_dims[k]) == gauss_binom(n,k,q) for k in range(d+1)]):
#             break
#     for i in range(d):
#         assert P_dim.count(i) == gauss_binom(n,i,q), "got {}, expected {}".format(P_dim.count(i),gauss_binom(n,i,q))
    relations = []
    P_list = list(P)
#     print(len(P))
#     print(P_list)
#     for k in range(d):
#         for U in P_by_dims[k]:
#             for W in P_by_dims[k+1]:
    dictionary = dict(zip(P_list, list(range(len(P_list)))))
    for i in range(d,0,-1):
#         print(i)
        for W in P_by_dims[i]:
            l = P_list.index(W)
            for V in W.subspaces(i-1):
                k = dictionary[V]
                relations.append([k, l])
#     for i, j in itertools.product(range(len(P)), repeat=2):
#         if (P_list[i].dimension() == (P_list[j].dimension() - 1)) and P_list[i].is_subspace(P_list[j]):
#             relations.append([i, j])
    elem_labs={}
    for W in P:
        elem_labs=[W.basis_matrix() for W in P_list]
    return Poset(data=(set(range(len(P))),relations)), elem_labs, F.base_ring().cardinality()

# An instantiation of the poset of linear subspaces of dimension <=d of an n-dimensional linear space over the field of q elements.
# Will include methods of a cell complex such as boundary and coboundary.

class grassmann_poset(object):
    def __init__(self, n, d, q):
        self.poset, self.elements_labels, char = build_grassmann_poset(n, d, q)
        self.field_size = q
        self.coboundary_coeficients = 3 if char == 2 else 2
        self.top_rank = d - 1

    #         for i in range(self.top_rank):
    #             assert not np.any(np.matmul(self.get_couboundary_matrix(i+1), self.get_couboundary_matrix(i)) % Gr.coboundary_coeficients),  "problem with the coboundary operators: delta_{}*delta_{} != 0".format(i+1,i)

    def show(self):
        self.poset.show(element_labels=self.get_elements_labels())

    def get_level_sets(self):
        return self.poset.level_sets()

    def get_level_sets_lengths(self):
        return [len(l) for l in self.poset.level_sets()]

    def get_element(self, i, j):
        return self.poset.level_sets()[i][j]

    def boundary(self, element):
        return self.poset.lower_covers(element)

    def coboundary(self, elements, i, f):
        target = self.poset.level_sets()[i + 1]
        l = len(target)
        g = np.zeros(l, dtype=int)
        for element in elements:
            S = self.poset.upper_covers(element)
            for i in range(l):
                if S.count(target[i]) > 0:
                    g[i] += 1
        return g % f

    def get_couboundary_matrix(self, i):
        domain = self.poset.level_sets()[i]
        target = self.poset.level_sets()[i + 1]
        d = len(domain)
        l = len(target)
        g = np.zeros(shape=(l, d), dtype=int)
        for k in range(d):
            S = self.poset.upper_covers(domain[k])
            for j in range(l):
                if S.count(target[j]) > 0:
                    g[j][k] += 1
        return g % self.coboundary_coeficients

    def get_hasse_diag(self):
        return self.poset.hasse_diagram()

    def get_elements_labels(self):
        return dict([(i, self.elements_labels[i]) for i in range(len(self.elements_labels))])

    def get_nx_graph(self):
        E = self.get_hasse_diag().edges()
        pos = self.get_hasse_diag().positions()
        return nx.relabel_nodes(nx.DiGraph([(e[0], e[1]) for e in E]), self.get_elements_labels())

    def to_xml(self, path):
        nx.write_graphml_lxml(self.get_nx_graph(), path)


# Example:
# Gr = grassmann_post(5,3,2)
# Gr.show()
