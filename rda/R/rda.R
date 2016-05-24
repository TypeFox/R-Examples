rda <- function(x, y, xnew=NULL, ynew=NULL, prior=table(y)/length(y), 
                alpha=seq(0, 0.99, len=10), delta=seq(0, 3, len=10),
                regularization="S", genelist=FALSE, trace=FALSE){
  
  this.call <- match.call()

  #if ( !((is.numeric(alpha))*(is.numeric(delta))*(is.vector(alpha))*
  #     (is.vector(delta))) ) {
  #  stop("Both parameters alpha and delta have to be numeric scalars or vectors.")
  #}

  #if ( !((is.numeric(x))*(is.numeric(y))*(is.matrix(x))*
  #       (is.vector(y))) ) {
  #  stop("Both x and y have to be numeric. See help for details.")
  #}   
  
  if ( regularization=="S" ){
    tmp <- rda1(x=x, y=y, xnew=xnew, ynew=ynew,
                prior=prior, alpha=alpha, delta=delta, 
                genelist=genelist, trace=trace)
  }
  else if ( regularization=="R" ){
    tmp <- Drda1(x=x, y=y, xnew=xnew, ynew=ynew,
                 prior=prior, alpha=alpha, delta=delta, 
                 genelist=genelist, trace=trace)
  }
  else {
    stop("Regularization argument must be either 'S' or 'R'. ")
  }

  tmp$call <- this.call; tmp$reg <- regularization
  class(tmp) <- "rda"
  tmp
}
