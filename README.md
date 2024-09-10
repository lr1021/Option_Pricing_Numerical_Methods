# Option_Pricing_Numerical_Methods

This repository contains implementations of numerical methods for option pricing, based on **Hull, J. C. (2003).** *Options, Futures, and Other Derivatives* (2nd ed.). Prentice-Hall, specifically Chapter 14: Numerical Procedures.

## Table of Contents
- [Overview](#overview)
- [Implemented Methods](#implemented-methods)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [References](#references)
- [License](#license)

## Overview

The purpose of this repository is to provide a practical implementation of various numerical methods for pricing options as described in Chapter 14: Numerical Procedures of John Hull's textbook *Options, Futures, and Other Derivatives*. These methods can be used to value derivative securities where analytical solutions may not be available or practical.
In addition to this, hedge parameters such as delta, gamma, theta, vega and rho can be estimated using prices obtained numerically (by recomputing the desired option price with a small change in relevant varaible) or with procedures specific to the used numerical method.

### Key Topics:
- Building Desired Option (To be implemented)
- Implicit Finite Difference
- Explicit Finite Difference (To be implemented)
- Binomial Trees
- Monte Carlo Simulation (To be implemented)
- Hedge Parameter Estimation (To be implemented)

## Implemented Methods

The following numerical methods have been implemented:

1. **Binomial Option Pricing Model (Cox-Ross-Rubinstein Model)**
   - A tree-based method to model the evolution of the underlying asset's price and calculate the option price.

2. **Trinomial Option Pricing Model**
   - An extension of the binomial model that provides more accuracy by considering three possible price movements (up, down, unchanged).

3. **Finite Difference Methods (FDM)**
   - Techniques to solve partial differential equations (PDEs) related to option pricing, including the explicit, implicit, and Crank-Nicolson methods.

4. **Monte Carlo Simulation**
   - A stochastic method used to simulate the underlying asset's price paths and calculate the option price by averaging the payoffs.

5. **Other Methods**
   - Any other methods included from Chapter 14, such as approximations and lattice models, will be documented here.

## Dependencies

To run the code in this repository, you need the following dependencies:

- Python 3.x
- NumPy
- Pandas
- Matplotlib (for plotting results)
- Jupyter (optional, for running notebooks)

You can install these dependencies using pip:

```bash
pip install numpy pandas matplotlib jupyter
```

## Usage

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/lr1021/Option_Pricing_Numerical_Methods.git
   cd Option_Pricing_Numerical_Methods
   ```

2. **Run the Notebooks or Scripts:**
   - You can explore the Jupyter notebooks provided in the `notebooks/` directory, which contain explanations, code, and examples. If there are rendering problems within the Jupyter notebooks, have a look at the homonymous .pdf files in the same directory
   - You can view the numerical method implementations in the Python scripts in the `methods/` directory, which contain class definitions and methods (with some helper functions mostly for visualisation)

3. **Customise Option to be priced and Method Parameters:**
   - Each method allows for customization of priced option and parameters like the number of time steps, volatility, risk-free rate, and more. Modify these parameters as needed within the scripts or notebooks.

## References

- **John C. Hull. (2003).** *Options, Futures, and Other Derivatives* (2nd ed.). Prentice-Hall.

This textbook is a widely recognized reference for derivatives pricing and includes detailed explanations of the numerical methods implemented in this repository.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

### Example Notebooks and Code:

Feel free to explore the provided notebooks that demonstrate how each method can be applied to price European and American options, as well as more complex derivatives.

---
