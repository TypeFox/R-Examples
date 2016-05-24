################################################
## WMTSA package generic functionality
##
##  Functions:
##
##    crystal.names
##    eda.plot
##    reconstruct
##    wavStackPlot
##
################################################

###
# crystal.names
##

"crystal.names" <- function(x, ...)
  UseMethod("crystal.names")

###
# eda.plot
##

"eda.plot" <- function(x, ...)
  UseMethod("eda.plot")

###
# reconstruct
###

"reconstruct" <- function(x,...)
  UseMethod("reconstruct")

###
# wavStackPlot
##

"wavStackPlot" <- function(x, ...)
  UseMethod("wavStackPlot")

  