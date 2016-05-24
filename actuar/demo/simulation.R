### ===== actuar: An R Package for Actuarial Science =====
###
### Demo of the portfolio simulation facilities provided by actuar
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

require(actuar)


## A simple Compound Poisson model: S_t = C_1 + ... + C_{N_t}, with
## N_t ~ Poisson(10), C ~ Lognormal(log(1500) - 1, 1). The names of
## the components serve no purpose here but are required.
pf <- simul(list(y = 10), model.freq = expression(y = rpois(10)),
            model.sev = expression(y = rlnorm(log(1500) - 1, 1)))
pf                                      # print method
aggregate(pf)                           # aggregate claim amounts
frequency(pf)                           # frequencies
severity(pf)                            # individual claim amounts
severity(pf, splitcol = 10)             # last period separate


## Simple (continuous) mixture of models: S_t|Theta ~ Poisson(Theta),
## Theta ~ Gamma(2, 1). Any names can be used in the model.
pf <- simul(list(Theta = 1, S = 10),
            model.freq = expression(Theta = rgamma(2, 1), S = rpois(Theta)))
aggregate(pf)                           # actual data
frequency(pf)                           # same, here

## Model with with mixtures for both frequency and severity.
pf <- simul(list(entity = 10, year = 5),
            model.freq = expression(entity = rgamma(2, 1),
                year = rpois(entity)),
            model.sev = expression(entity = rnorm(5, 1),
                year = rlnorm(entity, 1)))
pf
aggregate(pf)
frequency(pf)

## Same model as above, but with weights incorporated into the model.
## The string "weights" should appear in the model specification
## wherever weights are to be used.
wit <- runif(10, 2, 10)
(wit <- runif(50, rep(0.5 * wit, each = 5), rep(1.5 * wit, each = 5)))
(pf <- simul(list(entity = 10, year = 5),
             model.freq = expression(entity = rgamma(2, 1),
                 year = rpois(weights * entity)),
             model.sev = expression(entity = rnorm(5, 1),
                 year = rlnorm(entity, 1)),
             weights = wit))
weights(pf)                             # extraction of weights

## Three level hierarchical model (sector, unit, contract). Claim
## severity varies only by sector and unit. The number of "nodes" at
## each level is different.
nodes <- list(sector = 2, unit = c(3, 4),
              contract = c(10, 5, 8, 5, 7, 11, 4), year = 6)
mf <- expression(sector = rexp(2),
                 unit = rgamma(sector, 0.1),
                 contract = rgamma(unit, 1),
                 year = rpois(weights * contract))
ms <- expression(sector = rnorm(2, sqrt(0.1)),
                 unit = rnorm(sector, 1),
                 contract = NULL,
                 year = rlnorm(unit, 1))
wijkt <- runif(50, 2, 10)
wijkt <- runif(300, rep(0.5 * wijkt, each = 6), rep(1.5 * wijkt, each = 6))
pf <- simul(nodes, model.freq = mf, model.sev = ms, weights = wijkt)
frequency(pf)
weights(pf)
