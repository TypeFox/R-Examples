#########################################
## PHASE 1: Initiation.
##          Fast. Run on scheduler?
##
## NB: Number in file groups.txt used to determine Phase2 width.
##
## Two ways to call this:
##  0: use cross object as provided for perm test
##  1,2,...,nSim: run nSim times with each spawning multiple Phase2 runs.
##
## Input files:
##       cross.RData
##       params.txt
##
## Output files:
##       Phase1.RData
##       groups.txt

library(qtlhot)
parallel.qtlhot(1,commandArgs(TRUE))
