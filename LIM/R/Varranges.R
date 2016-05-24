
##==============================================================================
##  Calculates ranges of inverse variables
##==============================================================================

Varranges <- function(lim,...) {

  if (lim$NVariables == 0)
    return(NULL)
  res   <- varranges (lim$A, lim$B,
                      lim$G, lim$H,
                      lim$VarA, lim$VarB, ...)
  rownames(res)<-lim$Variables

  return(res)
}
