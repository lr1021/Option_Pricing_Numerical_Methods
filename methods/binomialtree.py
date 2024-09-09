import numpy as np
from .node import Node
from .option import Option


def to_subscript(number):
    # Define a dictionary for subscript characters
    subscript_map = {
        '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄',
        '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'}
    # Convert the number to a string and map each digit to its subscript
    return ''.join(subscript_map[digit] for digit in str(number))


class BinomialTree:
    def __init__(self, time, steps):
        # Time attributes
        self.time = time
        self.steps = steps
        self.timestep = self.time/self.steps
        # Following attributes need to be directly assigned or set with fit()
        # method.

        # Option to be evaluated by BinomialTree
        self.option = Option()
        # Interest rate r, continuous divided rate of asset q (for Futures
        # contracts set q = r)
        self.r = 0
        self.q = 0

        # Underlying asset price at time zero S, asset volatility sigma
        self.sigma = 0
        self.S = 1

        # Initialized value of f, the price of the option considered
        self.f = 0

        # Initialized up and down factors and up-probability of each Binomial
        # split
        self.u = 1
        self.d = 1/self.u
        self.p = 1

        # Initialized value of BinomialTree Nodes, stored in a dictionary
        self.nodes = {}
        for i in range(0, self.steps+1):
            for j in range(0, i+1):
                self.nodes[(i, j)] = Node(self.timestep*i,
                                          option=self.option,
                                          r=self.r,
                                          sigma=self.sigma)

    def __repr__(self):
        # Summary of Binomial Tree attributes
        return self.__class__.__name__ + "\n" + "(" + 'time ' + str(self.time)\
               + ' in ' + str(self.steps) + ' steps' + ")"\
               "\n" + "(" + 'r = ' + str(self.r) +\
               ', q = ' + str(self.q) + ', sigma = ' + str(self.sigma) + ")"\
               + '\n' + '(' + 'S = ' + str(self.S) + ', f = ' +\
               str(round(self.f, 4)) + ')'

    def __str__(self):
        # Tree representation at each level
        tree_str = []

        # Iterate over nodes
        for i in range(self.steps + 1):
            # Level representation
            level_str = []

            for j in range(i + 1):
                node = self.nodes.get((i, j))

                # Subscripts
                i_sub = to_subscript(i)
                j_sub = to_subscript(j)

                # Format the node as (S, f) and append to the level string list
                if node:
                    level_str.append(
                        f"({node.S:.2f}, {node.f:.2f}){i_sub}{j_sub}")
                else:
                    # Placeholder if empty
                    level_str.append("(None, None)")

            # Join the nodes in the level with spaces and center-align them
            tree_str.append(" ".join(level_str).center(20 * (self.steps + 1)))

            if i < self.steps:
                # Branches layer
                slashes = []
                for j in range(i + 1):
                    slashes.append(" /       \\")
                # Join the slashes and center-align them
                slashes_representation = "       ".join(slashes)\
                    .center(20 * (self.steps + 1))
                tree_str.append(slashes_representation)

        # Join levels with newlines and return the final tree representation
        return "\n".join(tree_str)

    def load_nodes(self):
        # Updates self.nodes with BinomialTree parameter values
        for i in range(0, self.steps+1):
            for j in range(0, i+1):
                self.nodes[(i, j)] = Node(self.timestep*i,
                                          option=self.option,
                                          r=self.r,
                                          sigma=self.sigma)

    def fit(self, S=None, r=None, q=None, sigma=None, option=None):
        """
        Calculates asset prices at all nodes in the Binomial tree using forward
        induction, then computes the option prices using backward induction.

        Parameters:
        ----------
        S : float
            Initial asset price at the root node.
        r : float
            Risk-free interest rate.
        q : float
            Asset dividend yield.
        sigma : float
            Volatility of the underlying asset.
        option : object
            An object representing the option, with methods to calculate
            the payoff and other properties.

        Steps:
        ------
        1. **Forward Induction:** Calculate asset prices at each node in the
             Binomial tree by iterating through node layers using up and down
             factors.

        2. **Backward Induction:** Compute the option prices by starting at the
             terminal nodes and moving backward to the root.
            - At terminal nodes, calculate the option payoff.
            - For each preceding node, calculate the option price using a
              discounted weighted average of the child nodes' prices.
            - If option.american is True then set option price the maximum
              between the above value and the option's intrinsic value (payoff
              at current node)

        Returns:
        --------
        None
            The function updates the asset prices and option prices at all
            nodes in the Binomial tree.
        """

        # if given new values, update BinomialTree attributes, nodes and
        # up/down/p factors
        if S:
            self.S = S
        if r:
            self.r = r
        if q:
            self.q = q
        if sigma:
            self.sigma = sigma
        if option:
            self.option = option
        self.load_nodes()

        self.u = np.e**(self.sigma*(self.timestep)**(1/2))
        self.d = 1/self.u
        self.p = (np.e**((self.r-self.q)*self.timestep)-self.d)/(self.u-self.d)

        # Forward induction to set asset prices in BinomialTree nodes
        for i in range(0, self.steps+1):
            for j in range(0, i+1):
                self.nodes[(i, j)].S = self.S*(self.u**j)*(self.d**(i-j))

        # Backward induction to set option prices in BinomialTree nodes
        for i in reversed(range(0, self.steps+1)):
            for j in range(0, i+1):
                self.nodes[(i, j)]
                # Terminal nodes
                if i == self.steps:
                    self.nodes[(i, j)].f = self.nodes[(i, j)].evaluate()
                    # print(f'{(i,j)}:{round(self.nodes[(i, j)].f,3)}')

                else:
                    # Discounted expected value
                    expected = (self.p*self.nodes[(i+1, j+1)].f +
                                (1-self.p)*self.nodes[(i+1, j)].f)\
                                    * np.e**(-self.r*self.timestep)
                    # American/European style option
                    if self.option.american:
                        self.nodes[(i, j)].f = max(
                            expected,
                            self.nodes[(i, j)].evaluate())
                    else:
                        self.nodes[(i, j)].f = expected

                    # print(f'{(i,j)}:{round(self.nodes[(i, j)].f,3)}')
        self.f = self.nodes[(0, 0)].f
