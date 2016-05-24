## ----installation, eval=FALSE--------------------------------------------
#  vignette("coala-install", package = "coala")

## ----create_model--------------------------------------------------------
library(coala)
model <- coal_model(sample_size = 3, loci_number = 1)

## ----print_model---------------------------------------------------------
model

## ------------------------------------------------------------------------
model <- model + feat_mutation(rate = 1, model = "IFS")
model

## ----echo=FALSE----------------------------------------------------------
funcs <- ls("package:coala")
funcs[grep("^feat_", funcs)]

## ------------------------------------------------------------------------
model <- coal_model(sample_size = c(5, 2, 0), loci_number = 1) +
  feat_migration(rate = 0.5, symmetric = TRUE) +
  feat_pop_merge(0.5, 3, 2) +
  feat_pop_merge(0.8, 2, 1)

## ------------------------------------------------------------------------
model <- coal_model(3, 1) +
  feat_mutation(rate = 1) +
  sumstat_seg_sites()
model

## ----echo=FALSE----------------------------------------------------------
funcs[grep("^sumstat_", funcs)]

## ------------------------------------------------------------------------
set.seed(123)
sumstats <- simulate(model)

## ------------------------------------------------------------------------
names(sumstats)

## ------------------------------------------------------------------------
sumstats$seg_sites[[1]]

## ------------------------------------------------------------------------
model <- model + locus_single(500)
model

## ------------------------------------------------------------------------
sumstats <- simulate(model)
sumstats$seg_sites[[1]]
sumstats$seg_sites[[2]]

## ------------------------------------------------------------------------
model <- model + locus_averaged(2, 750)
sumstats <- simulate(model)
length(sumstats$seg_sites)

## ------------------------------------------------------------------------
model <- coal_model(5, 1) +
  feat_mutation(rate = par_named("theta")) +
  sumstat_seg_sites()
sumstats1 <- simulate(model, pars = c(theta = 2.5))
sumstats2 <- simulate(model, pars = c(theta = 4.3))

## ----priors--------------------------------------------------------------
model <- coal_model(5, 1) +
  feat_mutation(rate = par_prior("theta", runif(1, 0, 10))) +
  sumstat_seg_sites()
sumstats <- simulate(model)
sumstats$pars
sumstats2 <- simulate(model)
sumstats2$pars

## ------------------------------------------------------------------------
model <- coal_model(5, 1) +
  feat_mutation(rate = par_range("theta", 0.1, 5)) +
  sumstat_seg_sites()

## ------------------------------------------------------------------------
model <- coal_model(4, 2) +
  feat_mutation(rate = par_named("theta")) +
  feat_recombination(rate = par_expr(theta * 2))

