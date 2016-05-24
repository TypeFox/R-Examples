## zzz.R (2013-11-27)

##   Library Loading

## Copyright 2013 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

.coalescentMCMCenv <- new.env()

assign("MCMCstats",
       data.frame(row.names = c("Number of trees output",
                  "Burn-in period", "Sampling frequency",
                  "Number of generations", "Nb of accepted moves")),
       envir = .coalescentMCMCenv)
