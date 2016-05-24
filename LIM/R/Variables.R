
##==============================================================================
##  Creates a list with the inverse variables
##==============================================================================

Variables <- function(lim, res=NULL, ...) {

  if (lim$NVariables == 0)
    return(NULL)

  if (is.null(res))
    res <- Lsei.lim(lim,parsimonious=TRUE,...)$X

  variables <- data.frame(values=lim$VarA%*%res-lim$VarB)

  rownames(variables)<-lim$Variables

  return(variables)

}
