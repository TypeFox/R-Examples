##
##  PURPOSE:   Pseudo goodness-of-fit test for a normal mixture
##             * generic function
##
##  AUTHOR:    Arnost Komarek
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   20/08/2009
##
##  FUNCTIONS: NMixPseudoGOF.R
##
## ==================================================================

## *************************************************************
## NMixPseudoGOF
## *************************************************************
NMixPseudoGOF <- function(x, ...)
{
  UseMethod("NMixPseudoGOF")
}  
