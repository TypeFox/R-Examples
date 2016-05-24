ComputeProbaPostBlock <- function(x, alpha, epsilon, delta) {
  beta <- alpha * delta + (1-alpha) * (1-delta)
  ord <- order(beta, decreasing = FALSE)
  beta <- beta[ord]
  # fb sont les valeurs des constantes pour chaque individus
  matSup <- matrix(0, ncol(x) , ncol(x)+1)
  matSup[upper.tri(matSup)] <- 1
  matInf <- 1-matSup
  # Proba per individuals
  partinf <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * delta[ord] 
  partsup <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * (1-delta[ord])   
  baseSup <- log(sweep(x[,ord], 2, partsup, "*") + sweep(1-x[,ord], 2, 1-partsup, "*"))
  baseInf <- log(sweep(x[,ord], 2, partinf, "*") + sweep(1-x[,ord], 2, 1 - partinf, "*"))
  fij <- exp(baseSup %*% matSup + baseInf%*% matInf)
  repere <- (c(beta, 1)-c(0, beta))
  probaInd <- fij %*% repere
  return(log(probaInd))
}



#' Computation of the model Cramer'v.
#' 
#' This function computes the model Cramer's V for a binary data set.
#' @param x a binary matrix.
#' @param param an instance of S4 class MvBinaryResult (provided by the function MvBinaryEstim)
#' @return Return the logprobability for each row of matrix x conditionally on the model defined by param.
#' 
#' @export
#' 
MvBinaryProbaPost <- function(x, param){
  if ( (is.matrix(x)==FALSE) || any((x==0) + (x==1) == 0) )
    stop("The input parameter x must be a binary matrix")
  if (class(param)!= "MvBinaryResult")
    stop("The imput parameter must be an instance of S4 class MvBinaryResult (provided by the function MvBinaryEstim) ")
  
  logproba <- rep(0, nrow(x))
  for (b in 1:max(param@blocks)){
    who <- which(param@blocks==b)
    if (length(who)>1){
      logproba <- logproba + ComputeProbaPostBlock(x[,who], param@alpha[who], param@epsilon[who], param@delta[who]) 
    }else{
      logproba <- logproba + x[,who]*log(param@alpha[who]) + (1-x[,who]) * log(1-param@alpha[who])
    }
  }
  return(logproba)
}