#########################################
## PHASE 5: Combine results for post-processing.
##          Fast. Run on scheduler.
##
## Input files (results of Phase 4):
##       Phase1.RData
##       mcmcN.RData (N = 1, ..., nruns)
##
## Output files:
##       Final.RData

library(qtlnet)
parallel.qtlnet(5)

