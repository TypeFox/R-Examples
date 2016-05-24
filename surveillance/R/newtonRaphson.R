################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Michaela's own implementation of a Newton-Raphson optimizer
###
### Copyright (C) 2010-2012 Michaela Paul
### $Revision: 589 $
### $Date: 2013-07-08 10:25:30 +0200 (Mon, 08 Jul 2013) $
################################################################################


#####################
# x - initial parameter values
# control arguments:
# scoreTol - convergence if max(abs(score)) < scoreTol
# paramTol - convergence if rel change in theta < paramTol
# F.inc - eigenvalues of the hessian are computed when the Cholesky factorization
#         fails, and a constant added to the diagonal to make the smallest
#         eigenvalue= F.inc * largest

# fn must return loglikelihood with score and fisher as attributes
#   fn <- function(theta,...){
#     ll <- loglik(theta,...)
#     attr(ll,"score") <- score(theta,...)
#     attr(ll,"fisher") <- fisher(theta,...)
#     return(ll)
#   }

newtonRaphson <- function(x,fn,..., control=list(), verbose=FALSE){

  # set default values
  control.default <- list(scoreTol=1e-5, paramTol=1e-8, F.inc=0.01,
                          stepFrac=0.5, niter=30)
  control <- modifyList(control.default, control)
  
  # number of step reductions, not positive definite Fisher matrices during iterations
  steph <- notpd <- 0
  convergence <- 99
  i <- 0
  
  rel.tol <- function(x,xnew){
    sqrt(sum((xnew-x)^2)/sum(x^2))
  }
  
  score <- function(fn){
    return(attr(fn,"score"))
  }
  fisher <- function(fn){
    return(attr(fn,"fisher"))
  }
  
  ll0 <- c(fn(x,...))
  if(verbose>1) cat("initial loglikelihood",ll0,"\n\n")
  
  # fn cannot be computed at initial par
  if(!is.finite(ll0) | is.na(ll0)){
    cat("fn can not be computed at initial parameter values.\n")
    return(list(convergence=30, notpd = notpd, steph = steph))
  }
      
  while(convergence != 0 & (i< control$niter)){
    i <- i+1
    ll <- fn(x,...)
    if(max(abs(score(ll))) < control$scoreTol){
      convergence <- 0
      break
    }
    # get cholesky decompositon
    F <- fisher(ll)
    F.chol <- try(chol(F),silent=TRUE)
    # could still give a nearly singular matrix
    # => could also check condition number
    
    if(inherits(F.chol,"try-error")){
      if(verbose>1) cat("fisher is not pd\n")
      # fisher is not pd
      notpd <- notpd +1
      ev <- eigen(F,symmetric=TRUE, only.values=TRUE)$values
      #  add a constant to diag(F)
      diag(F) <- diag(F) + (control$F.inc*(max(abs(ev))) - min(ev))/(1-control$F.inc)
      # compute cholesky decomposition of modified fisher
      F.chol <- chol(F)
    }
    direction <- chol2inv(F.chol)%*% score(ll)
    if(max(abs(direction)) < control$paramTol*(max(abs(x))+1e-8) ){
      convergence <- 0
      break
    }
    # do Newton-Raphson step
    x.new <- c(x + direction)
    ll.new <- fn(x.new,...)
    if(verbose>1) cat("iteration",i,"\trel.tol =",rel.tol(x,x.new),"\tabs.tol(score) =",max(abs(score(ll.new))),"\n")
    if(verbose>2) cat("theta =",round(x.new,2),"\n")
    if(verbose>1) cat("loglikelihood =",ll.new,"\n")
      

    ## Backtracking: reduce stepsize until we really improve the loglikelihood
    # ll(x1,lambda) = ll(x0) + lambda * fisher(x0)^-1 %*% score(x0)
    i.backstep <- 0
    ## Gray (2001) Ch 3: Unconstrained Optimization and Solving Nonlinear Equations
    # It is technically possible to construct sequences where ll(x1) > ll(x0)
    # at each step but where the sequence never converges.
    # For this reason a slightly stronger condition is usually used.
    # Dennis and Schnabel (1983): Numerical Methods for Unconstrained 
    #     Optimization and Nonlinear Equations. SIAM. (ch 6,3.2, p.126)
    # recommend requiring that lambda satisfy
    #   ll(x1) > ll(x0) + 1e-4 *(x1-x0)' %*% score(x0)
    while((is.na(ll.new) || (ll.new < c(ll)+ (1e-4)*sum(direction*score(ll)))) & (i.backstep <= 20)){
      if(verbose>1 & i.backstep==0) cat("backtracking: ")
      i.backstep <- i.backstep +1
      steph <- steph +1
      # reduce stepsize by a fixed fraction stepFrac
      direction <- control$stepFrac*direction
      x.new <- c(x + direction)
      ll.new <- fn(x.new,...)
      if(verbose>1) cat("*")
    }
    if(verbose & i.backstep>0) cat("\n")
    if(i.backstep >20){
      if(verbose>1)cat("backtracking did not improve fn\n")
      #cat("ll: ",ll,"\tll.new: ",ll.new,"\n")
      convergence <- 10
      break
    }

    x <- c(x.new)
      
    if(verbose>1) cat("\n")
  }
  ll <- fn(x,...)
  
  # max number of iterations reached, but check for convergence
  if(max(abs(score(ll))) < control$scoreTol){
    convergence <- 0
  }

  # convergence if
  # 1) relative difference between parameters is small
  # 2) absolute value of gradient is small
  # 3) stop after niter iterations
  
  if(i==control$niter & convergence !=0){
      if(verbose>1) cat("Newton-Raphson stopped after",i,"iterations!\n")
      # iteration limit reached without convergence
      convergence <- 10
  }
  if(verbose>1) cat("iteration",i,"\trel.tol =",rel.tol(x,x.new),"\tabs.tol(score) =",max(abs(score(ll))),"\n")
  if(verbose>2) cat("theta =",round(x.new,2),"\n")
  if(verbose>1) cat("loglikelihood =",c(ll),"\n\n")
 
  # loglikelihood
  loglik <- c(ll)
  
  # fisher info
  F <- fisher(ll)
  if(inherits(try(solve(F),silent=TRUE),"try-error")){ 
    cat("\n\n***************************************\nfisher not regular!\n")
    #print(summary(x))
    return(list(coefficients=x, loglikelihood=loglik, fisher=FALSE, convergence=22, notpd = notpd, steph = steph))
  }
  
  # check if solution is a maximum (i.e. if fisher is pd )
  eps <- 1e-10
  if(!all(eigen(F,symmetric=TRUE, only.values=TRUE)$values > eps)){
    if(verbose>1) cat("fisher information at solution is not pd\n")
    return(list(coefficients=x, loglikelihood=loglik, fisher=FALSE, convergence=21, notpd = notpd, steph = steph))
  }
 
  if(verbose>0) 
    cat("number of iterations = ",i," coverged = ", convergence ==0," log-likelihood = ",loglik, " notpd = ", notpd, " steph = ", steph, "\n")
  result <- list(coefficients=x, loglikelihood=loglik, fisher=FALSE, 
                 convergence=convergence, notpd=notpd, steph=steph,niter=i)
  return(result)
  
}
