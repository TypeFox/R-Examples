################ commands that have moved to RobExtremes

.mv2RobExtremes <- function(what){
  cat("\n\n", what, gettext(" has moved to package RobExtremes."),
      "\n\n")
}

kMAD <- function(x,k) .mv2RobExtremes("kMAD")
Pareto <- function(shape, Min) .mv2RobExtremes("Pareto")
GPareto <- function(loc, scale, shape, location) .mv2RobExtremes("GPareto")
GEV <- function(loc, scale, shape, location) .mv2RobExtremes("GEV")

EULERMASCHERONICONSTANT <- list(NULL)
class(EULERMASCHERONICONSTANT) <- "moved2RobExtremeClass"
APERYCONSTANT <- list(NULL)
class(APERYCONSTANT) <- "moved2RobExtremeClass"
print.moved2RobExtremeClass <- function(x,...)
   .mv2RobExtremes("Constants 'EULERMASCHERONICONSTANT' and 'APERYCONSTANT'")
