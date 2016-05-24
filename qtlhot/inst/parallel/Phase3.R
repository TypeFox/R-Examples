#########################################
## PHASE 3: Combine bic files into saved.scores.
##          Fast. Run on scheduler.
##
## Input files (results of Phase 1 and Phase 2):
##       Phase1.RData
##       permN.RData (N = 1, ..., n.split))
##
## Output files (one per invocation):
##       Phase3.RData
##
## See Phase 2 for explanation of bicN.RData files.

library(qtlhot)
parallel.qtlhot(3)
