#########################################
## PHASE 3: Combine bic files into saved.scores.
##          Fast. Run on scheduler.
##
## Input files (results of Phase 1 and Phase 2):
##       Phase1.RData
##       bicN.RData (N = 1, ..., nrow(groups))
##
## Output files (one per invocation):
##       Phase3.RData
##
## See Phase 2 for explanation of bicN.RData files.

library(qtlnet)
parallel.qtlnet(3)
