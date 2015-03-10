'''
A high-level class used for modeling various StateMatrix data-structures.
@author: Todd Gillette and Parsa Hosseini
'''

import math
import collections
import decimal

class AbstractMatrix():
    """
    A high-level matrix with capability of setting row and column names and
    dimensional information.
    """

    def __init__(self, nrows, ncols):
        """
        An abstract matrix capable of having row and column data points.
        """
        self.nrows = nrows
        self.ncols = ncols
        self.data = [[0.0 for _ in range(self.ncols)] 
                     for _ in range(self.nrows)]
        self.size = 0.0 # total sum of the matrix

    def get_width(self):
        """
        Retrieve how many columns there are in the current matrix.
        @return: matrix width (number of columns).
        """
        if len(self.data) == 0:
            return 0
        else:
            return len(self.data[0])
    
    def calculate_size(self):
        """
        Compute the total sum of the matrix by summing all row values.
        """
        self.size = float(sum([sum(col) for col in self.data]))

    def get_height(self):
        """
        Retrieve the number of rows comprising the current matrix.
        @return: matrix height (number of rows).
        """
        return len(self.data)

    def is_square(self):
        """
        Determines if all rows in the matrix are the same length.
        @return: boolean.
        """
        return len(set([len(row) for row in self.data])) == 1

class ContingencyMatrix(AbstractMatrix):
    """
    A Contingency Matrix (CM) is a model for quantifying categorical
    variables. Each Contingency Matrix is a 2x2 matrix, representing
    counts between two groups (A and B).
    """

    def __init__(self, node, i_g, i_not_g, not_i_g, not_i_not_g):
        """
        Creates a Contingency Matrix instances representing counts given
        counts within groups A and B.

        @param node: Domain of interest.
        """
        super(ContingencyMatrix, self).__init__(2, 2) # 2x2 matrix
        self.name = node
        # set values in the respective matrix
        self.data[0][0], self.data[0][1] = i_g, i_not_g
        self.data[1][0], self.data[1][1] = not_i_g, not_i_not_g
        self.calculate_size()

    def ipf(self):
        """
        Compute Iterative Proportional Fitting (IPF) on a contingency
        matrix in-place. Given a matrix with cells f(0, 0), f(1, 0), f(0, 1)
        and f(1, 1), cells f(0, 0) and f(1, 1) will have the same value, x.
        Values for f(1, 0) and f(0, 1) will however be (N/2) - x.
        """
        top = self.size * math.sqrt(self.get_g_x() * self.get_not_g_not_x())
        bottom = 2 * (math.sqrt(self.get_not_g_not_x() * self.get_g_x()) +
                      math.sqrt(self.get_not_g_x() * self.get_g_not_x()))
        self.data[1][1] = self.data[0][0] = top / bottom

        # Next, compute IPF for f(0,1) and f(1,0)
        diag_value = (self.size / 2) - self.get_g_x()
        self.data[1][0] = self.data[0][1] = diag_value
        self.calculate_size()

    ##########################################################################
    ###             Contingency Matrix value-access behaviors              ###
    ##########################################################################
    def get_g(self, to_prob=False):
        """
        Compute total probability or summation for group G (i.e. n(G)).

        @param to_prob: Output probability or raw count boolean.
        @return: Abundance of group G.
        """
        #val = float(self.__group_a.graph.total_count)
        val = float(self.data[0][0] + self.data[1][0])
        return (val / self.size) if to_prob else val

    def get_not_g(self, to_prob=False):
        """
        Compute total probability or summation for not-group G (i.e n(!G)).

        @param to_prob: Output probability or raw count boolean.
        @return: Abundance of not-group G.
        """
        #count = float(self.__group_b.graph.total_count)
        #return (count / self.size) if to_prob else count
        val = float(self.data[0][1] + self.data[1][1])
        return (val / self.size) if to_prob else val

    def get_g_x(self, to_prob=False):
        """
        Compute probability or count of variable X in group G (i.e. n(X, G).

        @param to_prob: Output probability or raw count boolean.
        @return: Abundance of X in group G.
        """
        val = float(self.data[0][0])
        return (val / self.size) if to_prob else val

    def get_not_g_not_x(self, to_prob=False):
        """
        Compute probability or count of variable not-X in not-group G
        (i.e. n(!X, G).

        @param to_prob: Output probability or raw count boolean.
        @return: Abundance of not-X in not-group G.
        """
        val = float(self.data[1][1])
        return (val / self.size) if to_prob else val

    def get_not_g_x(self, to_prob=False):
        """
        Compute probability or count of variable X in not-group G
        (i.e. n(X, !G).

        @param to_prob: Output probability or raw count boolean.
        @return: Abundance of X in not-group G.
        """
        val = float(self.data[0][1])
        return (val / self.size) if to_prob else val

    def get_g_not_x(self, to_prob=False):
        """
        Compute probability or count of variable not-X in group G
        (i.e. n(!X, G)).

        @param to_prob: Output probability or raw count boolean.
        @return: Numeric representing abundance of !X in group G.
        """
        val = float(self.data[1][0])
        return (val / self.size) if to_prob else val

    def get_x(self, to_prob=False):
        """
        Compute probability or count of variable X across both groups
        (i.e. n(X)).
        
        @param to_prob: Output probability or raw count boolean.
        @return: Abundance of X in groups G and !G.
        """
        val = float(self.data[0][0] + self.data[0][1])
        return (val / self.size) if to_prob else val

    def get_not_x(self, to_prob=False):
        """
        Compute probability or count of variable not-X across both groups
        (i.e. n(!X)).
        
        @param to_prob:  Output probability or raw count boolean.
        @return: Abundance of not-X in groups G and !G.
        """
        val = float(self.data[1][0] + self.data[1][1])
        return (val / self.size) if to_prob else val

    ##########################################################################
    ###                           Statistical Metrics                      ###
    ##########################################################################
    def get_lift(self):
        """
        Compute Lift (LI) metric.
        @return: Measure of over-representation
        """
        numerator = self.get_g_x(to_prob=True)  # P(AB)
        denominator = self.get_x(to_prob=True) * self.get_g(
            to_prob=True)  # {P(A) * P(B)}
        return numerator / denominator

    def get_cosine(self):
        """
        Compute Cosine (CO) metric.
        @return: Measure of over-representation
        """
        top = self.get_g_x(to_prob=True)
        bottom = math.sqrt(self.get_x(to_prob=True) * self.get_g(to_prob=True))
        return top / bottom

    def get_jaccard(self):
        """
        Compute Jaccard (JAC) metric.
        @return: Measure of over-representation.
        """
        return self.get_g_x(True) / (
            self.get_x(True) + self.get_g(True) - self.get_g_x(True))

    def get_laplace(self):
        """
        Compute Laplace (LP) correction.
        @return: Measure of over-representation.
        """
        return (self.get_g_x() + 1) / (self.get_x() + 2)

    def get_confidence(self):
        """
        Compute Confidence (CF) metric.
        @return: Measure of over-representation.
        """
        return self.get_g_x(to_prob=True) / self.get_x(to_prob=True)

    def get_phi_coefficient(self):
        """
        Compute Phi coefficient (PHI).
        @return: Measure of over-representation.
        """
        top = self.get_g_x(True) - (self.get_x(True) * self.get_g(True))
        bottom = math.sqrt(
            self.get_x(True) * self.get_g(True) * (1 - self.get_x(True)) * (
                1 - self.get_g(True)))
        return top / bottom

    def get_kappa_coefficient(self):
        """
        Compute Kappa coefficient (K).
        @return: Measure of over-representation.
        """
        top = self.get_g_x(True) + self.get_not_g_not_x(True) - self.get_x(
            True) * self.get_g(True) - self.get_not_x(True) * \
            self.get_not_g(True)
        bottom = 1 - self.get_x(True) * self.get_g(True) - self.get_not_x(
            True) * self.get_not_g(True)
        return top / bottom

    def get_hypergeometric_prob_mass(self):
        """
        Computes hypergeometric distribution given matrix.
        @return: p-value
        """
        k = self.get_g_x() # Number of target set domains that are the domain of interest
        N = self.size # Number of domains in target and baseline set
        n = self.get_g() # Number of domains in the target set
        K = self.get_x() # Occurances of the domain of interest in target and baseline sets
        n_give_n = combinatorial(N, n)

        m_give_x = combinatorial(K, k)
        nm_give_nk = combinatorial(N - K, n - k)
        dhyper = (m_give_x * nm_give_nk) / n_give_n
        return dhyper

    def get_hypergeometric_pval(self):
        """
        Computes hypergeometric distribution given matrix.
        @return: p-value
        """ 
        N = self.size # Number of domains in target and baseline set
        n = self.get_g() # Number of total domain occurences in the target set
        K = self.get_x() # Total occurences of the domain of interest (in target and baseline sets)
        k = self.get_g_x() # Number of target set domains that are the domain of interest
        n_give_n = combinatorial(N, n)
        phyper = 0

        max_range = min(K,n)
        if max_range - k > k:
            for x in range(0,int(k)):
                m_give_x = combinatorial(K, x)
                nm_give_nk = combinatorial(N - K, n - x)
                dhyper = (m_give_x * nm_give_nk) / n_give_n
                phyper += dhyper
            phyper = 1-phyper
        else:
            for x in range(int(k),int(max_range+1)):
                m_give_x = combinatorial(K, x)
                nm_give_nk = combinatorial(N - K, n - x)
                dhyper = (m_give_x * nm_give_nk) / n_give_n
                phyper += dhyper

        return phyper

    def measures(self):
        """
        For each Contingency Matrix, a pre-compiled dictionary of metrics
        will compute its TFBS abundance. This dictionary is such that the
        key and value is the metric and resultant output, respectively.

        :return: dictionary of metrics and its respective output.
        """
        measures = collections.OrderedDict()
        measures['LP'] = self.get_laplace()
        measures['CO'] = self.get_cosine()
        measures['JAC'] = self.get_jaccard()
        measures['LI'] = self.get_lift()
        measures['CF'] = self.get_confidence()
        measures['K'] = self.get_kappa_coefficient()
        measures['PHI'] = self.get_phi_coefficient()
        measures['p-value'] = self.get_hypergeometric_pval()
        measures['n(' + self.__group_a.name + '_baseline)'] = self.get_g_x()
        measures['n(' + self.__group_b.name + '_query)'] = \
            self.get_not_g_x()
        return measures

class StateMatrix(AbstractMatrix):
    ''' 
    A StateMatrix is your traditional multi-dimensional array whereby its 
    contents are indexed by respective integers.
    '''
    def __init__(self, nrows, ncols):
        super(StateMatrix, self).__init__(nrows, ncols)
        #self.state = [[0.0 for _ in range(self.ncols)] 
        #                        for _ in range(self.nrows)]
        self.transpose = MatrixTranspositionFactory(self)
        self.T = self.transpose

    def get_data(self,i,j):
        ''' 
        Get the data at a specific index.
        @param i: Row number
        @param j: Column number
        '''
        return self.data[i][j]
    
    def set_data(self,i,j,val):
        self.data[i][j] = val

    def transpose(self):
        '''
        Traditional transposition of a matrix.
        '''
        return self.transpose

    def debug(self):
        ''' 
        Useful method to print-out an entire matrix row-by-row.
        '''
        for rownum in range(self.nrows):
            print(self.data[rownum])
        #    print(self.data[rownum], self._state[rownum)]

class MatrixTranspositionFactory():
    '''
    Works with StateMatrix classes to provide low cost matrix transposition.
    '''
    
    def __init__(self,matrix):
        self.transpose = matrix
        self.T = matrix
        
    def get_data(self,i,j):
        return self.transpose.data[j][i]

    def set_data(self,i,j,val):
        self.transpose.data[j][i] = val

    def transpose(self):
        return self.T

class DirectionalMatrixWrapper():
    '''
    A wrapper for directional score matrices
    '''
    
    def __init__(self,nrows,ncols,T=None):
        self.nrows = nrows
        self.ncols = ncols
        if T is None:
            self.score = StateMatrix(nrows,ncols)
            self.extend_flag = StateMatrix(nrows,ncols)
            self.a_c_match = StateMatrix(nrows,ncols)
            self.T = DirectionalMatrixWrapper(nrows,ncols,self)
        else:
            self.T = T
            self.score = T.score.T
            self.extend_flag = T.extend_flag.T
            self.a_c_match = T.a_c_match.T

    def transpose(self):
        return self.T

def factorial(num):
    """
    Computes the factorial for a given integer. To enable factorial
    computation of large numbers, log-scaling summation is performed, with
    the resultant exponential of the value being returned.
    
    :param num: Input integer.
    :return: Factorial result.
    """
    log_scale = sum([math.log(i) for i in range(1, num + 1)])
    return decimal.Decimal(log_scale).exp()


def combinatorial(n, r):
    """
    Computes the combinatorial value given n-choose-r integers.
    :param n: Total size.
    :param r: Selection size.
    :return: combinatorial for n-choose-r items.
    """
    n, r = int(n), int(r) # counts are as floats; cast to integer
    numerator = factorial(n)
    denominator = factorial(r) * factorial(n - r)
    return numerator / denominator
