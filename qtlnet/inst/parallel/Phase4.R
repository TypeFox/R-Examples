#########################################
## PHASE 4: Sample Markov chain (MCMC).
##          Slow. Parallelize.
##
## NB: nruns is found in the params.txt parameter file.
##
## Input files (results of Phase 1 and Phase 2):
##       Phase1.RData
##       Phase3.RData
##
## Output files (one per invocation):
##       mcmcN.RData (N = 1, ..., nruns)
##
## Argument:
##       index (number from 1 to nruns)

library(qtlnet)
parallel.qtlnet(4, commandArgs(TRUE))
