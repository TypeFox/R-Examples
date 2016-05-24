#########################################
## PHASE 2: Compute BIC scores.
##          Slow. Parallelize.
##
## NB: Find number of lines in groups.txt file created in Phase 1.
##
## Input files (result of Phase 1):
##       Phase1.RData
##
## Output files: (one per invocation)
##       bicN.RData (N = 1, ..., nrow(groups)
##
## Argument:
##       index (1 to number of lines in groups.txt file)

library(qtlnet)
parallel.qtlnet(2, commandArgs(TRUE))
