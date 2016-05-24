#########################################################################
##
## hidden functions to help in model implementation
##
## NOTE: these are not exported to the user and should always be
##       used in model functions. As such, fixing problems here
##       fixes them in all functions simultaneously.
##
## ghislain.vieilledent@cirad.fr, March 2012
##
#########################################################################


##=======================================================================
##
## calling.function
## return name of the calling function
##
##=======================================================================

calling.function <- function(parentheses=TRUE) {
  calling.function <- strsplit(toString(sys.call(which=-3)),",")[[1]][1]
  if (parentheses){
    calling.function <- paste(calling.function, "()", sep="")
  }
  return(calling.function)
}

##=======================================================================
##
## Check group
##
##=======================================================================

check.group <- function(group,data) {

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

check.mcmc.parameters <- function(burnin, mcmc, thin) {
    
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
## Check interval
##
##=======================================================================

check.interval <- function (interval,nobs) {
  if (sum(!is.numeric(interval))!=0) {
    cat("Error: interval must be a numeric scalar or vector of length equal to the number of observations.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (length(interval)==1) {
    return(rep(interval,nobs))
  }
  if (length(interval)>1 & length(interval)!=nobs) {
    cat("Error: interval must be a numeric scalar or vector of length equal to the number of observations.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (length(interval)>1 & length(interval)==nobs) {
    return(interval)
  }
}

##=======================================================================
##
## Check verbose
##
##=======================================================================

check.verbose <- function (verbose) {
  if (!(verbose%in%c(0,1))) {
    cat("Error: verbose must take value 0 or 1.\n")
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

check.Y.Binomial <- function (Y) {
  if (sum(!(c(Y)%in%c(0,1)))>0) {
    cat("Error: Response variable must take value 0 or 1.\n")
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

form.seeds <- function(seed) {
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
## Check tuning parameters
##
##=======================================================================

check.tune <- function(tune) {
  if(length(tune)!=1 || !is.numeric(tune) || tune<=0) {
    cat("Error: tune must be a positive scalar.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(tune)
} 

##=======================================================================
##
## Check and form starting values
##
##=======================================================================

#===============================
# form beta.start
#===============================

form.beta.start <- function (fixed,data,beta.start,np,family,defaults=NA) {
 
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

form.sigma2.start <- function (fixed,data,sigma2.start,family) {
 
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

form.Vb.start <- function (Vb.start,nq) {

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

form.mvn.prior <- function(mubeta, Vbeta, np) {
  
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

form.wishart.prior <- function(r, R, nq) {
    
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

check.ig.prior <- function(nu, delta) {
     
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
