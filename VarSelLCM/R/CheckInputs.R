# Verifie les parametres d entrees
CheckInputs <- function(x, g, initModel, vbleSelec, discrim, paramEstim, nbcores, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep){
  if ( (is.numeric(g)==FALSE) || (length(g)!=1))
    stop("The component number have to be an integer of length one!")
  
  if ((is.data.frame(x)==FALSE) && (is.matrix(x)==FALSE))
    stop("Data set must be a data frame or a matrix!")
  
  if (is.logical(vbleSelec) == FALSE)
    stop("Input vbleSelec must be logical")
  
  if ((length(discrim) != ncol(x)) || (all(discrim %in% c(0,1))==FALSE))
    stop("Input discrim must be logical of length number of variables")
  
  if (is.logical(paramEstim) == FALSE)
    stop("Input paramEstim must be logicial")
  
  if (is.numeric(nbcores) == FALSE)
    stop("Input nbcores must be numeric")
  
  if ((is.numeric(nbSmall) == FALSE) || (length(nbSmall)!=1))
    stop("Input nbSmall must be numeric of size one")
  
  if ((is.numeric(iterSmall) == FALSE) || (length(iterSmall)!=1))
    stop("Input iterSmall must be numeric of size one")
  
  if ((is.numeric(nbKeep) == FALSE) || (length(nbKeep)!=1))
    stop("Input nbKeep must be numeric of size one")
  
  if ((is.numeric(iterKeep) == FALSE) || (length(iterKeep)!=1))
    stop("Input iterKeep must be numeric of size one")
  
  if ((is.numeric(tolKeep) == FALSE) || (length(tolKeep)!=1))
    stop("Input tolKeep must be numeric of size one")  
  

}