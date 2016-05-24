#########################################################################
##
## Following functions perform checks for MCMCpack hierarchical models
## (denoted MCMCh...).
## ghislain.vieilledent@cirad.fr, May 5 2011
##
#########################################################################

##=======================================================================
##
## Check group
##
##=======================================================================

check.group.hmodels <- function(group,data) {

  if (!(group %in% colnames(data))) {
    cat("Error: group must be a string which gives the name of the grouping variable in data \n")
    stop("Please respecify and call ", calling.function(), " again.",
    call.=FALSE)
  }
  return(0)
}

##=======================================================================
##
## Check mcmc parameters
##
##=======================================================================

check.mcmc.parameters.hmodels <- function(burnin, mcmc, thin) {
    
  if(mcmc %% thin != 0) {
    cat("Error: MCMC iterations not evenly divisible by thinning interval.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(mcmc <= 0) {
    cat("Error: MCMC iterations must be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE) 
  }
  if(burnin < 0) {
    cat("Error: Burnin iterations negative.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(((burnin+mcmc) %% 10 != 0) || (burnin+mcmc)<100) {
    cat("Error: Value 'burnin+mcmc' should be divisible by 10 and >= 100.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(thin < 1) {
    cat("Error: Thinning interval must be superior or equal to 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

##=======================================================================
##
## Check verbose
##
##=======================================================================

check.verbose.hmodels <- function (verbose) {
  if (!(verbose%in%c(0,1))) {
    cat("Error: verbose must take value 0 or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

##=======================================================================
##
## Check FixOD
##
##=======================================================================

check.FixOD.hmodels <- function (FixOD) {
  if (!(FixOD%in%c(0,1))) {
    cat("Error: FixOD must take value 0 or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

##=======================================================================
##
## Check Y in c(0,1) for Binomial process
##
##=======================================================================

check.Y.Binomial.hmodels <- function (Y) {
  if (sum(!(c(Y)%in%c(0,1)))>0) {
    cat("Error: Response variable must take value 0 or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

##=======================================================================
##
## Check Y is a positive integer for Poisson process
##
##=======================================================================

check.Y.Poisson.hmodels <- function (Y) {
  if (sum(!(c(Y)>=0 && c(Y)%%1==0)>0)) {
    cat("Error: Response variable must be a positive integer.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

##=======================================================================
##
## Check and form seed
##
##=======================================================================

form.seeds.hmodels <- function(seed) {
  if(length(seed)!=1) {
    cat("Error: Mersenne seed should be of length 1.\n")
    stop("Please respecify and call ", calling.function(), " again.", call.=FALSE)
  }
  if(length(seed)==1) {
    if(is.na(seed)) {
      seed <- 12345
    }
    else {
      seed <- as.integer(seed)
    }
    if(seed < 0) {
      cat("Error: Mersenne seed negative.\n")
      stop("Please respecify and call ", calling.function(), " again.", call.=FALSE)                       
    }
  }
  return(seed)
}

##=======================================================================
##
## Check and form starting values
##
##=======================================================================

#===============================
# form beta.start
#===============================

form.beta.start.hmodels <- function (fixed,data,beta.start,np,family,defaults=NA) {
 
  if (is.na(beta.start)[1] & is.na(defaults)[1]){ # use GLM estimates
    beta.start <- matrix(coef(glm(fixed, family=family, data=data)), np, 1)
  }
  else if(is.na(beta.start)[1] & !is.na(defaults)[1]){ # use passed values
    beta.start <- matrix(defaults,np,1)     
  }
  else if(is.null(dim(beta.start))) {
    beta.start <- beta.start * matrix(1,np,1)  
  }
  else if(!all(dim(beta.start) == c(np,1))) {
    cat("Error: beta.start not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
  }
  return(beta.start)
}

#===============================
# form sigma2.start
#===============================

form.sigma2.start.hmodels <- function (fixed,data,sigma2.start,family) {
 
  if (is.na(sigma2.start)[1]){ # use GLM estimates
    sigma2.start <- var(residuals(glm(fixed, family=family, data=data)))
  }
  else {
    sigma2.start <- as.integer(sigma2.start[1])
  }
  if (sigma2.start<=0) {
    cat("Error: Starting value for sigma2 negative.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  } 
  return(sigma2.start)
}

#===============================
# form Vb.start
#===============================

form.Vb.start.hmodels <- function (Vb.start,nq) {

  if (is.na(Vb.start)[1]) {
    Vb.start <- diag(1,nq)
  }
  else if(is.null(dim(Vb.start))) {
    Vb.start <- Vb.start * diag(nq)  
  }
  if ((dim(Vb.start)[1] != nq) || (dim(Vb.start)[2] != nq)) {
    cat("Error: Vb.start not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(diag(Vb.start)>0)!=nq) {
    cat("Error: Vb.start should have positive values on the diagonal.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(Vb.start)
}

##=======================================================================
##
## Check and form priors
##
##=======================================================================

#===============================
# form multivariate Normal prior
#===============================

form.mvn.prior.hmodels <- function(mubeta, Vbeta, np) {
  
  # prior mean
  if(is.null(dim(mubeta))) {
    mubeta <- mubeta * matrix(1,nrow=np,ncol=1)  
  } 
  if((dim(mubeta)[1] != np) || (dim(mubeta)[2] != 1)) {
    cat("Error: in N(mubeta,Vbeta) prior, mubeta not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
     
  # prior variance
  if(is.null(dim(Vbeta))) {
    if (length(Vbeta) > np){
      cat("Vbeta was passed as a vector longer than np.
           Vbeta must be either a scalar or a matrix.")
      stop("Please respecify and call ", calling.function(), " again.\n", call.=FALSE)
    }
    Vbeta <- Vbeta * diag(np)    
  }
  if((dim(Vbeta)[1] != np) || (dim(Vbeta)[2] != np)) {
    cat("Error: in N(mubeta,Vbeta), prior Vbeta not conformable [p times p].\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
  }

  # check Vbeta for symmetry
  symproblem <- FALSE
  for (i in 1:np){
    for (j in i:np){
      if (Vbeta[i,j] != Vbeta[j,i]){
        symproblem <- TRUE
      }
    }
  }
  if (symproblem){
    cat("Vbeta is not symmetric.\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)       
  }

  # check Vbeta for positive number on the diagonal
  if (sum(diag(Vbeta)>0)!=np) {
    cat("Error: Vbeta should have positive values on the diagonal.\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
  }   
  return(list(mubeta,Vbeta))
}

#===================
# form Wishart prior
#===================

form.wishart.prior.hmodels <- function(r, R, nq) {
    
  # check to see if degrees of freedom produces proper prior
  if(r < nq) {
    cat("Error: in Wishart(r,rR) prior, r less than q.\n")
    stop("Please respecify and call ", calling.function(), " again.\n")
  } 
    
  # form the prior scale matrix
  if(is.null(dim(R))) {
    R <- R * diag(nq)
  }
  if((dim(R)[1] != nq) | (dim(R)[2] != nq)) {
    cat("Error: in Wishart(r,rR) prior, R not comformable [q times q].\n")
    stop("Please respecify and call ", calling.function(), " again.\n")
  }

  # check R for symmetry
  symproblem <- FALSE
  for (i in 1:nq){
    for (j in i:nq){
      if (R[i,j] != R[j,i]){
        symproblem <- TRUE
      }
    }
  }
  if (symproblem){
    cat("R is not symmetric.\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)       
  }

  # check R for positive number on the diagonal
  if (sum(diag(R)>0)!=nq) {
    cat("Error: R should have positive values on the diagonal.\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
  }
  return(list(r,R))
}

#==========================
# check inverse Gamma prior
#==========================

check.ig.prior.hmodels <- function(nu, delta) {
     
  if(nu <= 0) {
    cat("Error: in IG(nu,delta) prior, nu less than or equal to zero.\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
  }
  if(delta <= 0) {
    cat("Error: in IG(nu,delta) prior, delta less than or equal to zero.\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)    
  }
  return(0)
}

#==========================
# END
#==========================
