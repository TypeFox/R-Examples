## ----sfs-----------------------------------------------------------------
sfs <- c(112, 57, 24, 34, 16, 29, 8, 10, 15)

## ----model setup---------------------------------------------------------
library(coala)
model <- coal_model(10, 50) +
  feat_mutation(par_prior("theta", runif(1, 1, 5))) +
  sumstat_sfs()

## ----simulate, cache=TRUE------------------------------------------------
sim_data <- simulate(model, nsim = 2000, seed = 17)

## ------------------------------------------------------------------------
# Getting the parameters
sim_param <- create_abc_param(sim_data, model)
head(sim_param, n = 3)

# Getting the summary statistics
sim_sumstat <- create_abc_sumstat(sim_data, model)
head(sim_sumstat, n = 3)

## ----abc, fig.align="center", fig.width=5--------------------------------
suppressPackageStartupMessages(library(abc))
posterior <- abc(sfs, sim_param, sim_sumstat, 0.05, method = "rejection")
hist(posterior, breaks = 20)

