## ------------------------------------------------------------------------
library(coala)
model <- coal_model(c(10, 15), 100) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_migration(par_range("m", 0, 3), symmetric = TRUE) +
  feat_pop_merge(par_range("t_split", 0.1, 2), 2, 1) + 
  feat_recombination(1) +
  sumstat_jsfs()

## ------------------------------------------------------------------------
library(jaatha)
jaatha_model <- create_jaatha_model(model)

