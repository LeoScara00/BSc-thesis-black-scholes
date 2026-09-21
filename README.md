# BSc Thesis: Black–Scholes Option Pricing

**Leonardo Scaramelli · Bocconi University**

This repository contains my BSc thesis on the Black–Scholes model and an R script for the revised numerical experiments.

## Files

- [BSc_Thesis.pdf](BSc_Thesis.pdf): the thesis, covering stochastic calculus, option pricing, and numerical simulation.
- [simulation.R](simulation.R): a self-contained script comparing exact simulation and Euler–Maruyama, studying Monte Carlo sampling convergence, and assessing time-discretisation error.

The experiments price a European call with initial stock price **24**, strike **26**, annual risk-free rate **3.29%**, annual volatility **25%**, and maturity **90/365 years**. The Black–Scholes benchmark is **0.538036793350**.

## Run the experiments

Only **base R** is required; no additional packages are needed. The verified run used R 4.1.2.

**In RStudio:** open `simulation.R`, select **Session → Set Working Directory → To Source File Location**, then click **Source**.

Alternatively, run this command in a terminal from the folder containing the script:

```sh
Rscript --vanilla simulation.R results
```

The script uses fixed random seeds and automatically creates a `results/` folder containing CSV results, LaTeX tables, PDF and PNG figures, and R session information. These outputs are generated locally and are not included in this repository. Re-running the script overwrites the generated outputs.
