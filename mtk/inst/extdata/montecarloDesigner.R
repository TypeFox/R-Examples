#' A basic Monte Carlo Designer
#' @param factors a vector of factors names
#' @param distribNames the name of a probability distribution or a vector of such names
#' @param distribParameters a list of lists of distribution parameters
#' @param size the required sample size
#' @param ...
#' @value a named list consisting of a design data.frame (main) and a list of information parameters (information)
#' @note this function uses the Quantiles function defined in the mtk package (file mtkCommonFunctions.R)
#' @examples
#'   MC <- basicMonteCarlo(LETTERS[1:4])
basicMonteCarlo <- function(factors, distribNames="unif", distribParameters=list(min=0,max=1), size, ...){

  ## PRELIMINARIES
  ## number of factors and their distributions
  nbf <- length(factors)
  if(length(distribNames) == 1){
    distribNames <- rep(distribNames, nbf)
    distribParameters <- rep(list(distribParameters), nbf)
  }
  
  ## MAIN CALCULATIONS
  N <- size*nbf
  design <- matrix(runif(N), nrow=size, ncol=nbf)
  ## quantile calculations
  for(i in seq(nbf)){
    design[,i] <- Quantiles(design[,i],	distribNames[i], distribParameters[[i]])
  }
  colnames(design) <- factors
  
  ## OUTPUT
  information <- list(SamplingMethod="Basic Monte Carlo")
  resultat <- list(main=as.data.frame(design), information=information)
  
  return(resultat)
}
