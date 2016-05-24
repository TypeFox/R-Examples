#########################################################################
##
## hidden functions to help in model implementation
##
## NOTE: these are not exported to the user and should always be
##       used in model functions. As such, fixing problems here
##       fixes them in all functions simultaneously.
##
## ghislain.vieilledent@cirad.fr, June 2011
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
    cat("Error: Burnin iterations must be positive.\n")
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

check.verbose <- function (verbose) {
  if (!(verbose %in% c(0,1))) {
    cat("Error: verbose must take value 0 or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

##=======================================================================
##
## Check save.rho
##
##=======================================================================

check.save.rho <- function (save.rho) {
  if (!(save.rho %in% c(0,1))) {
    cat("Error: save.rho must take value 0 or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

##=======================================================================
##
## Check save.p
##
##=======================================================================

check.save.p <- function (save.p) {
  if (!(save.p %in% c(0,1))) {
    cat("Error: save.p must take value 0 or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

##=======================================================================
##
## Check save.N
##
##=======================================================================

check.save.N <- function (save.N) {
  if (!(save.N %in% c(0,1))) {
    cat("Error: save.N must take value 0 or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

##=======================================================================
##
## Check data
##
##=======================================================================

check.T.poisson <- function (T,nobs) {
  if(length(T)!=nobs) {
    cat("Error: 'trials' must have the same length as the response variable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(!is.numeric(T)) {
    cat("Error: 'trials' must be a vector of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(T))>0) {
    cat("Error: 'trials' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(T!=0 & T!=1)>0) {
    cat("Error: 'trials' must be a vector of 0 or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

## check.Y.poisson.old <- function (Y,T) {
##   if(!is.numeric(Y)) {
##     cat("Error: 'counts' must be a vector of numeric values.\n")
##     stop("Please respecify and call ", calling.function(), " again.",
##          call.=FALSE)
##   }
##   if (sum(is.na(Y))>0) {
##     cat("Error: 'counts' must not contain missing values.\n")
##     stop("Please respecify and call ", calling.function(), " again.",
##          call.=FALSE)
##   }
##   if (sum(Y<0 | Y%%1!=0)>0) {
##     cat("Error: 'counts' must be a vector of positive integers.\n")
##     stop("Please respecify and call ", calling.function(), " again.",
##          call.=FALSE)
##   }
##   if (sum(Y!=0 & T==0)>0) {
##     cat("Error: 'counts' must be zero when 'trials' equals zero.\n")
##     stop("Please respecify and call ", calling.function(), " again.",
##          call.=FALSE)
##   }
##   return(0)
## }

check.Y.poisson <- function (Y) {
  if(!is.numeric(Y)) {
    cat("Error: 'counts' must be a vector of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(Y))>0) {
    cat("Error: 'counts' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(Y<0 | Y%%1!=0)>0) {
    cat("Error: 'counts' must be a vector of positive integers.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

check.T.binomial <- function (T,nobs) { #= This version of the function assumes T>0
  if(length(T)!=nobs) {
    cat("Error: 'trials' must have the same length as the response variable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(!is.numeric(T)) {
    cat("Error: 'trials' must be a vector of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(T))>0) {
    cat("Error: 'trials' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(T<=0 | T%%1!=0)>0) {
    cat("Error: 'trials' must be a vector of integers superior to zero.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

check.Y.binomial <- function (Y,T) {
  if(!is.numeric(Y)) {
    cat("Error: 'presences' must be a vector of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(Y))>0) {
    cat("Error: 'presences' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(Y<0 | Y%%1!=0)>0) {
    cat("Error: 'presences' must be a vector of positive integers.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  } 
  if (sum(Y>T)>0) {
    cat("Error: 'presences' must be less than or equal to 'trials'.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

check.U <- function (U,nobs) {
  if(length(U)!=nobs) {
    cat("Error: 'alteration' must have the same length as the response variable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(!is.numeric(U)) {
    cat("Error: 'alteration' must be a vector of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(U))>0) {
    cat("Error: 'alteration' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (!all(U>=0 & U<=1)) {
    cat("Error: 'alteration' must be in the interval [0,1].\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}


check.X <- function (X,n) {
  if(!is.numeric(c(X))) {
    cat("Error: 'suitability' only accept vectors of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(c(X)))>0) {
    cat("Error: 'suitability' do not accept vectors with missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (dim(X)[1]!=n) {
    cat("Error: Incorrect vector length for the 'suitability' argument.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

check.W <- function (W,n) {
  if(!is.numeric(c(W))) {
    cat("Error: 'observability' only accept vectors of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(c(W)))>0) {
    cat("Error: 'observability' do not accept vectors with missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (dim(W)[1]!=n) {
    cat("Error: Incorrect vector length for the 'observability' argument.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

check.neighbors <- function (n.neighbors,ncell,neighbors) {
  # Length of n.neighbors=ncell
  if(length(n.neighbors)!=ncell) {
    cat("Error: 'n.neighbors' must have a length equal to the number of cells.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  # Numeric values
  if(!is.numeric(neighbors)) {
    cat("Error: 'neighbors' must be a vector of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(!is.numeric(n.neighbors)) {
    cat("Error: 'n.neighbors' must be a vector of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  # NA
  if (sum(is.na(neighbors))>0) {
    cat("Error: 'neighbors' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(n.neighbors))>0) {
    cat("Error: 'n.neighbors' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  # Positive integer
  if (sum(neighbors<=0 | neighbors%%1!=0)>0) {
    cat("Error: 'neighbors' must be a vector of integers superior to zero.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(n.neighbors<0 | n.neighbors%%1!=0)>0) {
    cat("Error: 'n.neighbors' must be a vector of positive integers.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  # Number of neighbors inferior to ncell-1
  if (!all(n.neighbors < (ncell-1))) {
    cat("Error: 'n.neighbors' must not contain values superior to ncell-1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  # Number of neighbors and length of neighbors
  if (sum(n.neighbors)!=length(neighbors)) {
    cat("Error: 'neighbors' must be a vector of length equal to sum(n.neighbors).\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  # Check values in neighbors
  if (sum(!(neighbors %in% c(1:ncell)))>0) {
    cat("Error: 'neighbors' must be a vector of integers between 1 and ncell.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  # Check that the target cell is not in the list of neighbors: --> ToDoList
  return(0)
}

check.sites <- function (sites,nobs) {
  if(length(sites)!=nobs) {
    cat("Error: 'sites' must have the same length as the response variable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(!is.numeric(sites)) {
    cat("Error: 'sites' must be a vector of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(sites))>0) {
    cat("Error: 'sites' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(sites<=0 | sites%%1!=0)>0) {
    cat("Error: 'sites' must be a vector of integers superior to zero.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

check.cells <- function (cells,nsite) {
  if(length(cells)!=nsite) {
    cat("Error: 'spatial.entity' must be of length equals to the number of sites.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(!is.numeric(cells)) {
    cat("Error: 'spatial.entity' must be a vector of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(cells))>0) {
    cat("Error: 'spatial.entity' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(cells<=0 | cells%%1!=0)>0) {
    cat("Error: 'spatial.entity' must be a vector of integers superior to zero.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

##=======================================================================
##
## Check and form starting parameters for Metropolis-Hastings
##
##=======================================================================

form.beta.start <- function (beta.start,np) {
  if (sum(is.na(beta.start))>0) { 
    beta.start <- rep(0,np)
  }
  else if(!is.na(beta.start)[1] & length(beta.start)!=np) {
    beta.start <- rep(beta.start[1],np) 
  }
  else if(length(beta.start)!=np) {
    cat("Error: beta.start not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
  }
  return(beta.start)
}

form.gamma.start <- function (gamma.start,nq) {
  if (sum(is.na(gamma.start))>0) { 
    gamma.start <- rep(0,nq)
  }
  else if(!is.na(gamma.start)[1] & length(gamma.start)!=nq) {
    gamma.start <- rep(gamma.start[1],nq)  
  }
  else if(length(gamma.start)!=nq) {
    cat("Error: gamma.start not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
  }
  return(gamma.start)
}

check.Vrho.start <- function (Vrho.start) {
  if (length(Vrho.start)!=1) {
    cat("Error: Vrho.start should be a scalar.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (!is.numeric(Vrho.start)) {
    cat("Error: Vrho.start should be a numeric.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (Vrho.start<=0) {
    cat("Error: Vrho.start should be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(Vrho.start)
}

##=======================================================================
##
## Check and form priors
##
##=======================================================================

check.mubeta <- function(mubeta, np) {
  if (is.null(dim(mubeta))) {
    mubeta <- rep(mubeta,np) 
  }
  else if (length(mubeta)!=np) {
    cat("Error: mubeta not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(mubeta)
}

check.Vbeta <- function(Vbeta, np) {
  if (!all(Vbeta>0)) {
    cat("Error: Vbeta should be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (is.null(dim(Vbeta))) {
    Vbeta <- rep(Vbeta,np) 
  }
  else if (length(Vbeta)!=np) {
    cat("Error: Vbeta not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(Vbeta)
}

check.mugamma <- function(mugamma, nq) {
  if (is.null(dim(mugamma))) {
    mugamma <- rep(mugamma,nq) 
  }
  else if (length(mugamma)!=nq) {
    cat("Error: mugamma not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(mugamma)
}

check.Vgamma <- function(Vgamma, nq) {
  if (!all(Vgamma>0)) {
    cat("Error: Vgamma should be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (is.null(dim(Vgamma))) {
    Vgamma <- rep(Vgamma,nq) 
  }
  else if (length(Vgamma)!=nq) {
    cat("Error: Vgamma not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(Vgamma)
}

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

check.Vrho.max <- function (Vrho.max) {
  if (length(Vrho.max)!=1) {
    cat("Error: Vrho.max should be a scalar.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (!is.numeric(Vrho.max)) {
    cat("Error: Vrho.max should be a numeric.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (Vrho.max<=0) {
    cat("Error: Vrho.max should be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(Vrho.max)
}

form.priorVrho <- function (priorVrho) {
  if (is.numeric(priorVrho[1]) & priorVrho[1] > 0.0) {
    priorVrho <- priorVrho[1]
  }
  else if (priorVrho=="Uniform") {
    priorVrho <- -2.0
  }
  else if (priorVrho=="1/Gamma") {
    priorVrho <- -1.0
  }
  else {
    priorVrho <- -1.0
    cat("priorVrho has been set to \"1/Gamma\" \n")
  }
  return(priorVrho)
}

#==========================
# END
#==========================
