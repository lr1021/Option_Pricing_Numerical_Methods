

class Option:
    def __init__(self, american=False, payoff_vars=['S']):
        # These are the variables required by self.payoff, possible values are
        # in ['S', 'r', 'sigma', 'time']
        self.payoff_vars = payoff_vars

        # Payoff function
        self.payoff = lambda S: 1

        # american controls whether the option can be exercised early or not
        self.american = american


class Call(Option):
    def __init__(self, X, american=False, payoff_vars=['S']):
        # These are the variables required by self.payoff, possible values are
        # in ['S', 'r', 'sigma', 'time']
        self.payoff_vars = payoff_vars
        self.strike = X

        # Payoff function
        self.payoff = lambda S: max(0, S - self.strike)

        # american controls whether the option can be exercised early or not
        self.american = american


class Put(Option):
    def __init__(self, X, american=False, payoff_vars=['S']):
        # These are the variables required by self.payoff, possible values are
        # in ['S', 'r', 'sigma', 'time']
        self.payoff_vars = payoff_vars
        self.strike = X

        # Payoff function
        self.payoff = lambda S: max(0, self.strike - S)

        # american controls whether the option can be exercised early or not
        self.american = american
