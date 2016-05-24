#########################################
## PHASE 1: Initiation.
##          Fast. Run on scheduler.
##
## NB: File groups.txt used to determine Phase2 width.
##
## Input files:
##       cross.RData
##       params.txt
##
## Output files:
##       Phase1.RData
##       groups.txt

library(qtlnet)
parallel.qtlnet(1)
