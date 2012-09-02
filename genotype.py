"""
Contains the Genotype class.

"""

__author__ = 'Ricardo Azevedo'


import numpy as np
import numpy.random as rnd
import networkx as nx

class Genotype(object):
    """
    A gene network defined as described in: Wagner (1996) Does evolutionary plasticity evolve? *Evolution* 50: 1008-1023.

    **Attributes:**

    * **adjacency_matrix** -- adjacency matrix representation of the gene network
    * **connectivity** -- connectivity density
    * **n_genes** -- number of genes
    * **reg_matrix** -- matrix representation of the gene network such that columns represent regulators and rows represent their targets

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
        >>> net.adjacency_matrix
        array([[ 1.,  0.,  1.,  0.],
               [ 0.,  1.,  0.,  0.],
               [ 1.,  0.,  0.,  1.],
               [ 0.,  1.,  1.,  1.]])
        """
        assert type(matrix) is np.ndarray
        assert matrix.shape[0] == matrix.shape[1]
        self.reg_matrix = matrix

    @property
    def n_genes(self):
        return len(self.reg_matrix.diagonal())

    @property
    def adjacency_matrix(self):
        return (self.reg_matrix != 0) * 1.

    @property
    def connectivity(self):
        # multiplying by 1. changes the dtype from bool to float
        return self.adjacency_matrix.sum() / (self.n_genes * self.n_genes)

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
