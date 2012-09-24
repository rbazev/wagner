"""
genotype.py: implementation of the gene network model described in Siegal & Bergman (2002) Waddington's canalization
    revisited: developmental stability and evolution. PNAS 99: 10528-32.

Contains the Genotype class.

"""

import copy
import numpy as np
import numpy.random as rnd
#import pygraphviz as pgv
import networkx as nx
#from scipy.stats import poisson


__author__ = 'Ricardo Azevedo, Christina Burch, Kayla Peck, Amanda Whitlock'
__version__ = "0.0.1"

class Genotype(object):
    """
    A gene network defined as described in: Wagner (1996) Does evolutionary plasticity evolve? *Evolution* 50: 1008-1023.

    Attributes:

        adj_matrix: adjacency matrix representation of the gene network
        connectivity: connectivity density
        gene_network: matrix representation of the gene network such that columns represent regulators and rows represent
            their targets
        graph: PyGraphviz representation of the gene network
        n_genes: number of genes
        n_interactions: number of (nonzero) interactions in the gene network
        mean_abs_strength: mean absolute strength of interactions (excluding zeros)

    """

    def __init__(self, matrix):
        """
        Create Genotype object based on a matrix.  Check that the matrix is square and a numpy.ndarray.

        Parameters:

            matrix: matrix representation of the gene network such that columns represent regulators and rows represent
                their targets (type: numpy.ndarray)

        >>> net = Genotype(np.array([[-1.,  0., -1.,  0.], [ 0.,  3.,  0.,  0.], [ 1.,  0.,  0.,  1.], [ 0.,  2.,  1.,  1.]]))
        >>> net.connectivity
        0.5
        >>> net.n_genes
        4
        """
        assert type(matrix) is np.ndarray
        assert matrix.shape[0] == matrix.shape[1]
        self.gene_network = matrix

    @property
    def n_genes(self):
        """Number of genes in the network."""
        return len(self.gene_network.diagonal())

    @property
    def graph(self):
        """NetworkX representation of the gene network."""
        g = nx.DiGraph(data = self.gene_network.transpose())
        return g

    @property
    def interactions(self):
        """List of interactions in the form (regulator, target)."""
        return self.graph.edges()

    @property
    def n_interactions(self):
        """Number of nonzero interactions between genes."""
        return len(self.interactions)

    @property
    def connectivity(self):
        """Connectivity density of the network."""
        return  self.n_interactions / float(self.n_genes * self.n_genes)

#    @property
#    def mean_abs_strength(self):
#        """Mean absolute strength of interactions (excluding zeros)."""
#        return np.abs(self.gene_network).sum() / self.n_interactions
#
#    @property
#    def graph(self):
#        """PyGraphviz representation of the gene network"""
#        g = pgv.AGraph(directed = True, strict = False)
#        g.node_attr['fontname'] = 'helvetica'
#        g.node_attr['fontsize'] = 16
#        g.node_attr['shape'] = 'circle'
#        g.node_attr['color'] = 'gray'
#        g.node_attr['style'] = 'filled'
#        # add genes
#        for gene in range(self.n_genes):
#            g.add_node(gene)
#        # add interactions
#        for regulator in range(self.n_genes):
#            for target in range(self.n_genes):
#                weight = self.gene_network[target, regulator]
#                if weight:
#                    g.add_edge(str(regulator), str(target))
#                    e = g.get_edge(str(regulator), str(target))
#                    e.attr['penwidth'] = np.abs(weight) / self.mean_abs_strength
#                    e.attr['color'] = 'red'
#                    if weight < 0:
#                        e.attr['color'] = 'blue'
#                        e.attr['arrowhead'] = 'tee'
#        g.layout(prog = 'dot')
#        return g
#
#    def draw_graph(self, filename):
#        """Draw gene network using graphviz.  Output PNG file."""
#        self.graph.draw(filename)

    @property
    def connected_components(self):
        """
        Connected components.  Return list containing lists of nodes in each connected component.
        """
        return nx.connected_components(self.graph.to_undirected())

    @property
    def is_connected(self):
        """Whether the gene network is connected (type: bool)."""
        return len(self.connected_components) == 1

    def draw_graph(self):
        """
        Draw gene network using matplotlib.

        Note:
            These are crude plots.  Genes are labelled from 0.  Thick ends represent arrows.
            Self interactions are not shown.
        """
        nx.draw(self.graph)

    @staticmethod
    def generate_random(n_genes, connectivity):
        """
        Generate random genotype with a given n_genes and connectivity.  Sample weights from standard normal distribution.
        Ensure that connectivity is between 0 and 1.

        Note:
            The realized connectivity will be different from the entered value if n_genes * n_genes * connectivity is not
            an integer (see code example below).

        Parameters:
            n_genes: number of genes (type: int)
            connectivity: connectivity density (type: float)

        >>> net = Genotype.generate_random(4, .3)
        >>> net.connectivity
        0.25
        """
        assert 0 <= connectivity <= 1
        n_sites = n_genes * n_genes
        n_nonzero_sites = n_sites * connectivity
        flat_matrix = np.zeros(n_sites)
        flat_matrix[0:(n_nonzero_sites)] = rnd.normal(size = n_nonzero_sites)
        rnd.shuffle(flat_matrix)
        return Genotype(np.reshape(flat_matrix, (n_genes, n_genes)))

    def mutate(self, regulator, target, new_strength):
        """
        Replaces interaction strength by new value.

        Parameters:

            regulator: number of regulator gene
            target: number of target gene
            new_strength: new interaction strength
        """
        self.gene_network[target, regulator] = new_strength

    def mutate_random(self, n_mutations):
        """
        Mutate a given number of random interactions.  Sample new interaction strengths from standard normal distribution.
        Ensure that (i) mutations can only affect nonzero interactions, and (ii) the same interaction cannot be changed
        twice.

        Parameters:

            n_mutations: number of mutations
        """
        if n_mutations > self.n_interactions:
            n_mutations = self.n_interactions
        new_strengths = rnd.normal(size = n_mutations)
        mut_interactions = []
        while len(mut_interactions) < n_mutations:
            x = rnd.randint(0, self.n_interactions)
            if not x in mut_interactions:
                mut_interactions.append(x)
        for i in range(n_mutations):
            interaction = self.interactions[mut_interactions[i]]
            self.mutate(interaction[0], interaction[1], new_strengths[i])

    def set_mutation_rate(self, mutation_rate):
        """
        Set mutation rate per genome per generation, that is, the expected number of mutations per generation.
        """
        if mutation_rate > self.n_interactions:
            self.mutation_rate = self.n_interactions
        else:
            self.mutation_rate = mutation_rate

    def generate_asexual_offspring(self):
        """
        Generate copy of a Genotype, allowing mutations to occur.
        Return a new Genotype object.
        """
        offspring = copy.deepcopy(self)
        offspring.mutate_random(rnd.poisson(offspring.mutation_rate))
        return offspring

#    @staticmethod #not sure if staticmethod is correct for this
#    def mutate_genotype(n_genes, connectivity, mut_rate, matrix):
#        """
#        Mutates an existing matrix given n_genes, connectivity, mut_rate, and the gene network matrix (e.g. x.gene_network)
#
#        Note:
#            I am not sure if all of these need to be arguments or whether we can get the values from the earlier properties in the code
#            This method currently changes the matrix that is sent in to it - it does not keep the original matrix intact and create another one.
#            We'll have to discuss which way would be the best for our purposes.
#
#            Right now, the code reads through each element in the matrix and if it is non-zero, it mutates via a Poisson process
#            with probability mut_rate/(connectivity*n_sites). I need to check to make sure this is correct. If a mutation occurs,
#            the matrix element is replaced with an independent standard normal random variate.
#
#            An alternative way of writing this code would be to sum the number of non-zero elements, determine how many will mutate
#            using a Poisson process, then randomly choose the non-zero element(s) to mutate. I think this method would take longer
#            since you would need to loop through the matrix twice: once for counting non-zero elements and again to determine which
#            ones will mutate. But if there is some benefit to that method, it would be easy to switch.
#
#            CB: In making the decision between these alternatives, we should consider making the locations of non-zero matrix elements
#            an attribute of the genotype.
#
#        Parameters:
#            n_genes: number of genes (type: int)
#            connectivity: connectivity density (type: float)
#            mut_rate: mutation rate per individual network per generation (type: float)
#            matrix: gene network matrix
#
#        >>> net = Genotype.generate_random(4, .3)
#        >>> net.connectivity
#        0.25
#        >>> net_mut = Genotype.mutate_genotype(net.n_genes, net.connectivity, .1, net.gene_network)
#        >>> net_mut.connectivity
#        0.25
#
#        """
#        n_sites = n_genes * n_genes
#        for i in range(0,n_genes):
#            for j in range(0,n_genes):
#                if matrix[i][j] != 0:
#                    mutate = poisson.rvs(mut_rate/(connectivity*n_sites), size=1)
#                    #CB: Discuss whether this needs to be a random uniform variate between 0 and 1...
#                    if mutate > 0:
#                        #CB: and then this line would change to if mutate > mut_rate/(connectivity*n_sites):
#                        matrix[i][j] = rnd.normal(size=1)
#                        #print 'mutated element [' + repr(i) + '][' + repr(j) + ']'
#                    else:
#                        pass
#                else:
#                    pass
#        return Genotype(matrix)  #caution: changes original matrix!

if __name__ == "__main__":
    import doctest
    doctest.testmod()
