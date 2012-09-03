"""
Contains the Genotype class.

"""

__author__ = 'Ricardo Azevedo'

import numpy as np
import numpy.random as rnd
import pygraphviz as pgv
import networkx as nx

class Genotype(object):
    """
    A gene network defined as described in: Wagner (1996) Does evolutionary plasticity evolve? *Evolution* 50: 1008-1023.

    **Attributes:**

    * **adj_matrix** -- adjacency matrix representation of the gene network
    * **connectivity** -- connectivity density
    * **gene_network** -- matrix representation of the gene network such that columns represent regulators and rows represent their targets
    * **graph** -- PyGraphviz representation of the gene network
    * **n_genes** -- number of genes
    * **n_interactions** -- number of (nonzero) interactions in the gene network
    * **mean_abs_strength** -- mean absolute strength of interactions (excluding zeros)

    """

    def __init__(self, matrix):
        """
        Create Genotype object based on a matrix.  Check that the matrix is square and a *numpy.ndarray*.

        :param matrix: matrix representation of the gene network such that columns represent regulators and rows represent their targets
        :type matrix: numpy.ndarray
        :return: Genotype object

        >>> net = Genotype(np.array([[-1.,  0., -1.,  0.], [ 0.,  3.,  0.,  0.], [ 1.,  0.,  0.,  1.], [ 0.,  2.,  1.,  1.]]))
        >>> net.connectivity
        0.5
        >>> net.n_genes
        4
        >>> net.adj_matrix
        array([[ 1.,  0.,  1.,  0.],
               [ 0.,  1.,  0.,  0.],
               [ 1.,  0.,  0.,  1.],
               [ 0.,  1.,  1.,  1.]])
        """
        assert type(matrix) is np.ndarray
        assert matrix.shape[0] == matrix.shape[1]
        self.gene_network = matrix

    @property
    def n_genes(self):
        return len(self.gene_network.diagonal())

    @property
    def adj_matrix(self):
        """A matrix with the same structure as the **gene_network** but where 1 and 0 correspond to the presence and absence of an interaction, respectively."""
        # multiplying by 1. changes the dtype from bool to float
        return (self.gene_network != 0) * 1.

    @property
    def connectivity(self):
        """The connectivity density of the network."""
        return self.adj_matrix.sum() / (self.n_genes * self.n_genes)

    @property
    def n_interactions(self):
        return self.connectivity * np.power(self.n_genes, 2)

    @property
    def mean_abs_strength(self):
        """Mean absolute strength of interactions (excluding zeros)."""
        return np.abs(self.gene_network).sum() / self.n_interactions

    @property
    def graph(self):
        """PyGraphviz representation of the gene network"""
        g = pgv.AGraph(directed = True, strict = False)
        g.node_attr['fontname'] = 'helvetica'
        g.node_attr['fontsize'] = 16
        g.node_attr['shape'] = 'circle'
        g.node_attr['color'] = 'gray'
        g.node_attr['style'] = 'filled'
        # add genes
        for gene in range(self.n_genes):
            g.add_node(gene)
        # add interactions
        for regulator in range(self.n_genes):
            for target in range(self.n_genes):
                weight = self.gene_network[target, regulator]
                if weight:
                    g.add_edge(str(regulator), str(target))
                    e = g.get_edge(str(regulator), str(target))
                    e.attr['penwidth'] = np.abs(weight) / self.mean_abs_strength
                    e.attr['color'] = 'red'
                    if weight < 0:
                        e.attr['color'] = 'blue'
                        e.attr['arrowhead'] = 'tee'
        g.layout(prog = 'dot')
        return g

    def draw_graph(self, filename):
        """Draw gene network using graphviz.  Output PNG file."""
        self.graph.draw(filename)

    def connected(self):
        """
        Test whether gene network is connected.

        """
        ug = nx.Graph(data = self.graph)
        cc = nx.connected_components(ug)
        return len(cc) == 1

    @staticmethod
    def generate_random(n_genes, connectivity):
        """
        Generate random genotype with a given **n_genes** and **connectivity**.  Sample weights from standard normal distribution.

        .. note::
           The realized **connectivity** will be different from the entered value if **n_genes** * **n_genes** * **connectivity** is not an integer (see code example below).

        :param n_genes: number of genes
        :type n_genes: int
        :param connectivity: connectivity density
        :type connectivity: float
        :return: Genotype object

        >>> net = Genotype.generate_random(4, .3)
        >>> net.connectivity
        0.25
        """
        n_sites = n_genes * n_genes
        n_nonzero_sites = n_sites * connectivity
        flat_matrix = np.zeros(n_sites)
        flat_matrix[0:(n_nonzero_sites)] = rnd.normal(size = n_nonzero_sites)
        rnd.shuffle(flat_matrix)
        return Genotype(np.reshape(flat_matrix, (n_genes, n_genes)))



if __name__ == "__main__":
    import doctest
    doctest.testmod()
