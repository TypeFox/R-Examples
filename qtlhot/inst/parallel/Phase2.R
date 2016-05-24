#########################################
## PHASE 2: Compute BIC scores.
##          Slow. Parallelize.
##
## NB: Run this the number of times indicated in groups.txt from Phase1.
##
## Input files (result of Phase 1):
##       Phase1.RData
##
## Output files: (one per invocation)
##       permN.RData (N = 1, ..., n.split)
##
## Argument:
##       index (1 to number in groups.txt file)

library(qtlhot)
parallel.qtlhot(2, commandArgs(TRUE))
