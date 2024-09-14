from .node import Node
from .option import Option
import numpy as np
import pandas as pd

default_index_color = '#d9e1ff'


# Function that will be used in ImplicitFiniteDifference.show() to visualize
# model
def display(df, index_color=default_index_color, heat=False, scale='linear'):
    # Rename axes
    df = df.rename_axis('Asset Price (dollars)', axis='index') \
        .rename_axis('Time (years)', axis='columns')

    # Formatting and styling
    df_styled = df.style \
        .set_properties(**{'background-color': 'white'}) \
        .format("{:.2f}")

    def style_index(index):
        return f"background-color: {index_color};"

    df_styled = df_styled.map_index(style_index, axis=0)

    if heat:
        if scale == 'log':
            # Apply a logarithmic scale (handle zero values)
            gmap = np.log1p(df)  # log1p is log(1 + x) to handle zero values
        elif scale == 'quadratic':
            # Apply a quadratic scale
            gmap = np.square(df)
        else:
            # Default is linear scale
            gmap = df

        df_styled = df_styled.background_gradient(cmap='coolwarm', gmap=gmap,
                                                  axis=None)

    return df_styled


class ImplicitFiniteDifference:
    def __init__(self, time, tsteps, Smax, vsteps, fmax=0, fmin=0):

        # Time attributes
        self.time, self.tsteps = time, tsteps
        self.timestep = self.time/self.tsteps

        # Following attributes need to be directly assigned or set with fit()
        # method.

        # Option to be evaluated by ImplicitFiniteDifference
        self.option = Option()

        # Interest rate r, continuous divided rate of asset q (for Futures
        # contracts set q = r)
        self.r, self.q = 0, 0

        # Smax is a boundary value of S for which option price is known and
        # constant
        self.Smax = Smax
        self.vsteps = vsteps

        # Underlying asset price at time zero S and asset volatility sigma
        self.S = Smax/vsteps
        self.sigma = 0

        # valstep is the stock price difference considered in the model;
        # note this does not necessarily give self.S as a stock price
        # evaluated by model, .set_Smax() method is then used to adjust Smax
        # such that self.S is evaluated
        self.valstep = Smax/self.vsteps

        # fmax and fmin are option price for boundary conditions of the
        # underlying asset price S
        self.fmax = fmax
        self.fmin = fmin

        # Initialized value of f, the price of the option considered
        self.f = 0

        # This array will be used for computation through implicit finite
        # difference equation
        self.map = np.zeros((self.vsteps+1, self.tsteps+1))
        self.map[0, :] = self.fmax
        self.map[self.vsteps, :] = self.fmin

        # Initialized value of Nodes, stored in a dictionary
        self.nodes = {}
        for i in range(0, self.tsteps+1):
            for j in range(0, self.vsteps+1):
                self.nodes[(i, j)] = Node(self.timestep*i,
                                          option=self.option,
                                          r=self.r,
                                          sigma=self.sigma)

    def price_index(self, price_precision=2):
        # To be used in .show() method
        index = [f'{j*self.valstep:.{price_precision}f} \n [{j}]'
                 for j in reversed(range(0, self.vsteps+1))]
        return index

    def time_index(self, time_precision=2):
        # To be used in .show() method
        index = [f'{i*self.timestep:.{time_precision}f} \n [{i}]'
                 for i in range(0, self.tsteps+1)]
        return index

    def show(self, index_color=default_index_color, price_precision=2,
             time_precision=2, heat=False, scale='linear'):

        index = self.price_index(price_precision)
        columns = self.time_index(time_precision)

        map = pd.DataFrame(self.map, index=index, columns=columns)
        return display(map, index_color, heat=heat, scale=scale)

    def det_Smax(self):
        for k in range(1, self.vsteps):
            valstep = self.S/k
            # print(f'valstep = {valstep}')

            if self.Smax == self.S + (self.vsteps - k)*valstep:
                # Perfect fit
                self.valstep = valstep
                return

            elif self.Smax < self.S + (self.vsteps - k)*valstep:
                # valstep is too large
                pass

            elif self.Smax > self.S + (self.vsteps - k)*valstep:
                if k == 1:
                    # With the largest possible value of valstep, Smax cannot
                    # be evaluated by model
                    raise ValueError(f"Value of Smax {self.Smax} is too large,\
                                     decrease Smax or increase vsteps")
                else:
                    self.valstep = self.S/(k-1)
                    self.Smax = self.S + (self.vsteps - (k-1))*self.S/(k-1)
                    # valstep * vsteps larger smaller than Smax,
                    # use previous value and set self.Smax to smallest multiple
                    # of valstep  greater than Smax
                    return

    def load_map(self):
        # Update map with given values of fmax and fmin
        self.map = np.zeros((self.vsteps+1, self.tsteps+1))
        self.map[0, :] = self.fmax
        self.map[self.vsteps, :] = self.fmin

    def load_nodes(self):
        # Updates self.nodes with ImplicitFiniteDifference parameter values
        for i in range(0, self.tsteps+1):
            for j in range(0, self.vsteps+1):
                self.nodes[(i, j)] = Node(self.timestep*i,
                                          option=self.option,
                                          r=self.r,
                                          sigma=self.sigma)
                self.nodes[(i, j)].S = self.valstep*j
                self.nodes[(i, j)].evaluate()

    def show_intrinsic(self, reload=True, index_color=default_index_color,
                       price_precision=2, time_precision=2, heat=False,
                       scale='linear'):
        if reload:
            self.load_nodes()
        intrinsic = np.zeros((self.vsteps+1, self.tsteps+1))
        for i in range(0, self.tsteps+1):
            for j in range(0, self.vsteps+1):
                intrinsic[self.vsteps - j, i] = self.nodes[(i, j)].intrinsic

        index = self.price_index(price_precision)
        columns = self.time_index(time_precision)

        intrinsic_map = pd.DataFrame(intrinsic, index=index, columns=columns)

        return display(intrinsic_map, index_color, heat=heat, scale=scale)

    def fit(self, S=None, Smax=None, fmax=None, fmin=None, r=None, q=None,
            sigma=None, option=None):

        # if given new values, update ImplicitFiniteDifference attributes
        if S:
            self.S = S
        if Smax:
            self.Smax = Smax
            self.det_Smax()
        if fmax:
            self.fmax = fmax
        if fmin:
            self.fmin = fmin
        if r:
            self.r = r
        if q:
            self.q = q
        if sigma:
            self.sigma = sigma
        if option:
            self.option = option

        self.det_Smax()
        self.load_map()
        self.load_nodes()

        a = [1/2*self.r*j*self.timestep - 1/2*self.sigma**2*j**2*self.timestep
             for j in range(0, self.vsteps+1)]
        b = [1 + self.r*self.timestep + self.sigma**2*j**2*self.timestep
             for j in range(0, self.vsteps+1)]
        c = [-1/2*self.r*j*self.timestep - 1/2*self.sigma**2*j**2*self.timestep
             for j in range(0, self.vsteps+1)]

        for i in reversed(range(0, self.tsteps+1)):
            map = self.map.copy()
            if i == self.tsteps:
                for j in range(1, self.vsteps):
                    payoff = self.nodes[(i, j)].intrinsic
                    self.map[self.vsteps - j, i] = payoff
            else:
                fnext = map[1:-1, i+1]
                fnext[0] += -c[self.vsteps-1]*self.fmax
                fnext[-1] += -a[1]*self.fmin

                M = np.zeros((self.vsteps-1, self.vsteps-1))
                for row in range(0, self.vsteps-1):
                    if row > 0:
                        M[row, row-1] = c[self.vsteps-row-1]
                    M[row, row] = b[self.vsteps-row-1]
                    if row < self.vsteps-2:
                        M[row, row+1] = a[self.vsteps-row-1]
                Minv = np.linalg.inv(M)
                fcurr = Minv @ fnext

                if self.option.american:
                    for j in range(1, self.vsteps):
                        payoff = self.nodes[(i, j)].intrinsic
                        self.map[self.vsteps - j, i] =\
                            max(payoff, fcurr[self.vsteps-j-1])
                else:
                    self.map[1:self.vsteps-1, i] = fcurr

        self.f = self.map[self.vsteps-int(self.S/self.valstep), 0]

    def explicit_fit(self, S=None, Smax=None, fmax=None, fmin=None, r=None,
                     q=None, sigma=None, option=None):

        # if given new values, update ImplicitFiniteDifference attributes
        if S:
            self.S = S
        if Smax:
            self.Smax = Smax
            self.det_Smax()
        if fmax:
            self.fmax = fmax
        if fmin:
            self.fmin = fmin
        if r:
            self.r = r
        if q:
            self.q = q
        if sigma:
            self.sigma = sigma
        if option:
            self.option = option

        self.det_Smax()
        self.load_map()
        self.load_nodes()

        a = [((j*self.timestep)/(1+self.r*self.timestep)) *
             (-self.r+j*self.sigma**2)/2
             for j in range(0, self.vsteps+1)]
        b = [(1/(1+self.r*self.timestep))*(1-self.sigma**2*j**2*self.timestep)
             for j in range(0, self.vsteps+1)]
        c = [((j*self.timestep)/(1+self.r*self.timestep)) *
             (self.r+j*self.sigma**2)/2
             for j in range(0, self.vsteps+1)]

        for i in reversed(range(0, self.tsteps+1)):
            map = self.map.copy()
            if i == self.tsteps:
                for j in range(1, self.vsteps):
                    payoff = self.nodes[(i, j)].intrinsic
                    self.map[self.vsteps - j, i] = payoff
            else:
                fnext = map[1:-1, i+1]
                # fnext[0] += -c[self.vsteps-1]*self.fmax
                # fnext[-1] += -a[1]*self.fmin

                M = np.zeros((self.vsteps-1, self.vsteps-1))
                for row in range(0, self.vsteps-1):
                    if row > 0:
                        M[row, row-1] = c[self.vsteps-row-1]
                    M[row, row] = b[self.vsteps-row-1]
                    if row < self.vsteps-2:
                        M[row, row+1] = a[self.vsteps-row-1]
                fcurr = M @ fnext
                fcurr[0] += c[self.vsteps-1]*self.fmax
                fcurr[-1] += a[1]*self.fmin

                if self.option.american:
                    for j in range(1, self.vsteps):
                        payoff = self.nodes[(i, j)].intrinsic
                        if payoff > 0:
                            self.map[self.vsteps - j, i] =\
                                max(payoff, fcurr[self.vsteps-j-1])
                        else:
                            self.map[self.vsteps - j, i] =\
                                fcurr[self.vsteps-j-1]
                else:
                    self.map[1:self.vsteps-1, i] = fcurr

        self.f = self.map[self.vsteps-int(self.S/self.valstep), 0]
