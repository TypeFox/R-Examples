## JAGS data file for the normal model

## Observation count in phase II
k.II <- 10

## Mean and precision parameters for the priors
theta.mu <- c(0, 2, 0, 0)
theta.tau <- c(1, 1, 8, 8)
eta.mu <- c(0, 2, 0, 0)
eta.tau <- c(1, 1, 8, 8)

## Sample size of each observation
n.II <- c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10)

## Doses
d.II <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

## Observed responses, taken from the efficacy and safety models using the parameter values theta = eta = c(0, 2, 1, 1) (rounded to two decimals).
YE.II <- c(0.18, 0.33, 0.46, 0.57, 0.67, 0.75, 0.82, 0.89, 0.95, 1.00)
YS.II <- c(0.18, 0.33, 0.46, 0.57, 0.67, 0.75, 0.82, 0.89, 0.95, 1.00)

## Standard deviations
sigmaE <- 1
sigmaS <- 1

## Number of phase III trials
k.III <- 2
