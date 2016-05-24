##
##  PURPOSE:   Drawing traceplots of selected parameters from fitted objects
##             * generic function
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   16/01/2010
##
##  FUNCTIONS: tracePlots (16/01/2010)
##
## ======================================================================

## *************************************************************
## tracePlots
## *************************************************************
tracePlots <- function(x, ...)
{
  UseMethod("tracePlots")
}  
