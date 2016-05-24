###########################################################
# object of class character (= estimation method)
IRT.irfprob.character <- function(object, A, B, 
                                  xsi, theta, 
                                  guess = NULL, iIndex = 1:dim(A)[1],
                                  nnodes = nrow(theta), maxK = dim(A)[2],
                                  AXsi = matrix(0, nrow = dim(A)[1], ncol = maxK), ...){
  
  if(object %in% c("tam.mml", "tam.mml.2pl",
                   "tam.mml.mfr")){		
    res <- calc_prob.v5(iIndex, A, AXsi, B, xsi, theta, nnodes, 
                 maxK, recalc = TRUE)$rprobs
    
  }
  if(object %in% c("tam.mml.3pl")){	
    res <- .mml.3pl.calc_prob.v5(iIndex, A, AXsi, B, xsi, theta, 
                               nnodes, maxK, recalc = TRUE, guess)$rprobs
  }
  attr(res,"theta") <- theta
  return(res)
}
###########################################################


###########################################################
# object of class tam (tam.mml)
IRT.irfprob.tam <- function( object , ... ){
  ll <- object$rprobs
  dimnames(ll)[[1]] <- colnames(object$resp)
  attr(ll,"theta") <- object$theta
  attr(ll,"prob.theta") <- object$pi.k
  attr(ll,"G") <- object$G
  return(ll)
}
IRT.irfprob.tam.mml <- IRT.irfprob.tam
# IRT.irfprob.tam.mfr <- IRT.irfprob.tam		
###########################################################

###########################################################
# object of class tam (tam.mml)
IRT.irfprob.tam.mml.3pl <- function( object , ... ){
  ll <- object$rprobs
  dimnames(ll)[[1]] <- colnames(object$resp)
  attr(ll,"theta") <- object$theta
  attr(ll,"prob.theta") <- object$pi.k
  res <- list( "delta" = object$delta , 
               "delta.designmatrix" = object$delta.designmatrix )
  attr(ll,"skillspace") <- res	
  attr(ll,"G") <- object$G
  return(ll)
}
###########################################################

###########################################################
# objects of class tamaan
IRT.irfprob.tamaan <- function( object , ... ){
  if (object$tamaanify$method %in% c( "tam.mml" , "tam.mml.2pl")  ){
    res0 <- IRT.irfprob.tam( object , ... )			
  }
  if (object$tamaanify$method == "tam.mml.3pl"){
    res0 <- IRT.irfprob.tam.mml.3pl( object , ... )			
  }
  return(res0)	
}
###################################################################			