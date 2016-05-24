glmlep <-
function(x, y, family = c("gaussian", "binomial"), lambda=NULL,
           lambda.min=ifelse(n<p,0.05,0.001), nlambda=100, lambda2=0,
           kappa=ifelse(n<p,0.1,0.05), pen.fac=rep(1,p), tol=1e-6, max.ite = 1e3){
  
  family = match.arg(family)
  this.call = match.call()

  penalty.factor = pen.fac
  
  n <- dim(x)[1]
  p <- dim(x)[2] 
  if (length(y) != n) 
    stop(paste("number of observations in y (", length(y), ") not equal to the number of rows of x (", 
               n, ")", sep = ""))
  
  ## generate lambda sequence
  if(!is.null(lambda)){
    if (any(lambda < 0)) 
      stop("lambdas should be non-negative")
    nlambda <- length(lambda)
    ls <-  lambda
  }  else {
    if (lambda.min >= 1) 
      stop("lambda.min should be less than 1")
    ls <- SetLambda(x, y, lambda.min, nlambda, penalty.factor)
  }
  
  beta0<-rep(0,p+1)
  ## give estimates for different lambda values.
  beta = switch(family, 
                gaussian = sapply(1:length(ls), function(i) .C("gaulep",as.double(x),as.double(y),as.double(ls[i]),
                                                               as.double(kappa),as.double(tol),as.integer(max.ite),
                                                               as.integer(n),as.integer(p),as.double(beta0),DUP=FALSE)[[9]]),
                binomial = sapply(1:length(ls), function(i) .C("binlep",as.double(x),as.double(y),as.double(ls[i]),
                                                               as.double(kappa),as.double(tol),as.integer(max.ite),
                                                               as.integer(n),as.integer(p),as.double(beta0),DUP=FALSE)[[9]]))
  
  
  ## choose the one giving the minimum BIC value.
  loss <- sapply(1:nlambda, function(i)loglike(x, y, beta=beta[,i], family))
  df <- apply(beta[-1,], 2, function(x) sum(x!=0)) 
  
  EBIC <- loss + (log(log(p))+log(n))*df
  ind_l <- which.min(EBIC) 
  beta_min <- beta[,ind_l]
  
  fit <- list(beta=beta, lambda=ls, df=df, loss=loss, EBIC=EBIC, lambda.min=ls[ind_l], beta.min=beta_min)
  fit$call <- this.call
  class(fit) = c(paste(substr(family,1,3),"lep",sep=""), "glmlep")
  return(fit)
}
