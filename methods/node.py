from .option import Option


class Node:
    # Node objects make up the BinaryTree

    def __init__(self, time, S=1, option=Option(), r=0, sigma=0):
        self.S = S
        self.f = 0
        self.intrinsic = 0
        # Variables needed to compute payoffs are stored in the node.
        # These are can be defined uniquely for each node meaning interest rate
        # and volatility can be defined as varying over time (for example as
        # a function of time)
        self.time = time
        self.option = option
        self.r = r
        self.sigma = sigma

    def evaluate(self):
        # Variables are evaluated with node values
        variables = {var: getattr(self, var) for var in
                     self.option.payoff_vars}

        # Payoff is computed
        payoff = self.option.payoff(**variables)
        self.intrinsic = payoff
        return payoff
