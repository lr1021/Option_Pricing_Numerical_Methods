import numpy as np
import pandas as pd


def fit_cov(cov, vol, corr):
    if isinstance(cov, np.ndarray):
        numvars = cov.shape[0]
        vol = np.zeros(numvars)
        corr = np.zeros((numvars, numvars))
        for i in range(numvars):
            for j in range(numvars):
                vol[j] = cov[j, j]
                corr[i, j] = cov[i, j]/(vol[i]*vol[j])

    elif isinstance(corr, np.ndarray) and isinstance(vol, np.ndarray):
        numvars = len(vol)
        cov = np.zeros((numvars, numvars))
        for i in range(numvars):
            for j in range(numvars):
                covij = corr[i, j] * vol[i] * vol[j]
                cov[i, j] = covij

    elif isinstance(vol, np.ndarray):
        numvars = len(vol)
        corr = np.eye(numvars)
        cov = np.zeros((numvars, numvars))
        for i in range(numvars):
            cov[i, i] = vol[i]**2
    else:
        raise ValueError('Must define either cov, corr and vol, or vol' +
                         'as np.ndarrays')
    return cov, vol, corr


class run:
    def __init__(self, time, tsteps, starting, drift, vol):
        self.time = time
        self.tsteps = tsteps
        self.timestep = time/tsteps
        self.starting = starting
        self.drift = drift
        self.vol = vol

        self.values = 0
        self.differences = 0

    def simulate(self):
        values = np.zeros(self.tsteps+1)
        values[0] = self.starting
        differences = np.zeros(self.tsteps)
        random_elements = np.random.randn(self.tsteps)

        for step in range(self.tsteps):
            eps = random_elements[step]
            last_value = values[step]
            difference = last_value*(np.exp((self.drift - self.vol**2/2) *
                                            self.timestep +
                                            self.vol*self.timestep**(1/2) *
                                            eps) - 1)
            differences[step] = difference
            values[step+1] = last_value + difference
        self.values = values
        self.differences = differences


class multirun:
    def __init__(self, time, tsteps, starting, drift, cov='not_given',
                 vol='not_given', corr='not_given', var_names=None):
        self.time = time
        self.tsteps = tsteps
        self.timestep = time/tsteps

        numvars = len(starting)
        if var_names:
            self.var_names = var_names
        else:
            self.var_names = [str(i) for i in range(numvars)]

        self.starting = starting
        self.drift = drift
        self.cov, self.vol, self.corr = fit_cov(cov, vol, corr)

        self.values = 0
        self.differences = 0

    def simulate(self):
        numvars = len(self.starting)
        values = np.zeros((self.tsteps+1, numvars))
        values[0, :] = self.starting
        differences = np.zeros((self.tsteps, numvars))
        random_elements = np.random.multivariate_normal(
            np.zeros(numvars), self.cov, size=self.tsteps)

        for step in range(self.tsteps):
            eps = random_elements[step, :]
            last_values = values[step, :]
            step_difference = np.zeros(numvars)
            for var in range(numvars):
                vol = np.sqrt(self.cov[var, var])
                drift = self.drift[var]
                difference = last_values[var]*(np.exp((drift - vol**2/2) *
                                                      self.timestep +
                                                      self.timestep**(1/2) *
                                                      eps[var]) - 1)
                differences[step, var] = difference
                step_difference[var] = difference

            values[step+1, :] = last_values + step_difference
        self.values = values
        self.differences = differences
        self.final_values = values[-1, :]


class MonteCarlo:
    def __init__(self, time, tsteps, starting, drift, cov='not_given',
                 vol='not_given', corr='not_given',
                 var_names=None, payoff=lambda x: 1, nruns=100, r=0):

        self.time = time
        self.tsteps = tsteps
        self.timestep = time/tsteps
        self.nruns = nruns

        self.starting = starting
        self.drift = drift
        self.r = r

        numvars = len(starting)
        if var_names:
            self.var_names = var_names
        else:
            self.var_names = [str(i) for i in range(numvars)]
        self.cov, self.vol, self.corr = fit_cov(cov, vol, corr)

    def run(self, nruns=None):
        if nruns:
            self.nruns = nruns

        simulations = []
        model = multirun(self.time, self.tsteps, self.starting, self.drift,
                         cov=self.cov, vol=self.vol, corr=self.corr,
                         var_names=self.var_names)

        for run in range(self.nruns):
            model.simulate()
            simulations.append(pd.DataFrame(model.values,
                                            columns=self.var_names))
        self.simulations = simulations

    def price(self, payoff_func, variables):
        all_discounted_payoffs = []
        for simulation in self.simulations:
            if 'r' in self.var_names:
                mean_r = np.mean(simulation['r'])
            else:
                mean_r = self.r

            # Variables are evaluated with node values
            variables = {var: list(simulation[var])[-1] for var in variables}

            # Payoff is computed
            payoff = payoff_func(**variables)
            discounted_payoff = payoff * np.exp(-mean_r*self.time)
            all_discounted_payoffs.append(discounted_payoff)
        price = np.mean(all_discounted_payoffs)
        return price
