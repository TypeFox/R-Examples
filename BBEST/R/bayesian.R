#####################################################################################
#
#
#   FUNCTIONS THAT IMPLEMENT VARIOUS BAESIAN MODELS
#
#

#####################################################################################
# get.posterior(par, data, skel, knots.x, Gr, p.bkg)
#
# returns:      wrapper to the target function.
#               called from DEoptim
# arguments
#   par:        DEoptim parameter
#   data:       an object of type data
#   skel:       DEoptim parameter
#   knots.x:    the spline knot x-values.
#   Gr:         low-r G(r) information
#   p.bkg:      the probability that a single pixel contains "only" background.
get.posterior <- function(par, data, skel=NA, knots.x, Gr=NA, p.bkg=.5){
  par <- relist(par, skel)
 
  knots.y <- if(is.null(par$knots.y)) NA else par$knots.y
  alpha <- if(is.null(par$alpha)) 1 else par$alpha
  ADP <- if(is.null(par$ADP)) NA else par$ADP
  knots.n <- length(knots.x)
  pars <- if(is.null(par$pars)) NA else par$pars

  if(!is.na(ADP)){
#   recalculate coherent baseline
    n.atoms <- data$fitADP$n.atoms
	ADP <- rep(ADP, length(n.atoms))
	scatter.length <- data$fitADP$scatter.length
    N_total <- sum(n.atoms)
    f.av2 <- (sum(n.atoms*scatter.length)/N_total)^2
    f2.av <- sum(n.atoms*scatter.length^2)/N_total
	L <- (f.av2-f2.av)/f.av2
    expADP <- 0
    for(i in 1:length(n.atoms))
      expADP <- expADP + n.atoms[i]*scatter.length[i]^2*exp(-ADP[i]*data$x^2)/N_total/f2.av
    
	data$SB <- 1-expADP*(1-L)   
  }
  
  if(!is.na(pars[1]))
    posterior <- logPosteriorAnalyt(data=data, pars=pars, p.bkg=p.bkg) 
  else{
    posterior <- logPosterior(data=data, alpha=alpha, knots.x=knots.x, 
                              knots.y=knots.y, Gr=Gr, p.bkg=p.bkg)
  }
  return(posterior)
} 



####################################################################################
# logPosterior(data, alpha, knots.x, knots.y, Gr, p.bkg)
#
# returns:      the optimization function from Fischer et al.: the negative log of
#               the posterior.
# arguments
#   data:       an object of type data. See reference manual.
#   alpha:      scale parameter.
#   knots.x:    the spline knot x-values.
#   knots.y:    the spline knot y-values ('c' in Fischer et al.).
#   Gr:         low-r G(r) information
#   p.bkg:      the probability that a single pixel contains "only" background.
logPosterior <- function(data, alpha, knots.x, knots.y, Gr=NA, p.bkg=0.5) {
  Phi <- basisMatrix(x=data$x, knots.x=knots.x)	  
  # 1. Prior 
  psi.prior <- logPriorF(knots.x, knots.y)$prior

  # 2. Likelihood:
  psi.likelihood.bkg    <- logLikelihoodBkg(   knots.y=knots.y, y=data$y-data$SB, Phi=Phi, p.bkg=p.bkg, sigma=data$sigma)
  psi.likelihood.signal <- logLikelihoodSignal(knots.y=knots.y, y=data$y-data$SB, Phi=Phi, p.bkg=p.bkg, sigma=data$sigma, lambda=data$lambda)

  max.log <- pmax(psi.likelihood.bkg, psi.likelihood.signal)  # parallel max --> vector of max values
  #  avoiding inf values:
  psi.likelihood <- max.log + log(exp(psi.likelihood.bkg - max.log) + exp(psi.likelihood.signal - max.log))
  psi.likelihood <- -sum(psi.likelihood)

  # 3. Gr=-4*Pi*rho*r restriction
  psi.gr <- 0
  if(!is.na(Gr[1]))
    psi.gr <- logProbabilityBkgR(y=data$y-data$SB, Gr=Gr, alpha=alpha, Phi=Phi, knots.y=knots.y)
  
  # 4. Together
  psi <- psi.prior + psi.likelihood + psi.gr 
  return(psi)
}


logPosteriorAnalyt <- function(data, pars, p.bkg=0.5) {

  # 2. Likelihood:
  psi.likelihood.bkg <- logLikelihoodBkgAnalyt(pars=pars, x=data$x, y=data$y-data$SB, p.bkg=p.bkg, sigma=data$sigma)
  psi.likelihood.signal <- logLikelihoodSignalAnalyt(pars=pars, x=data$x, y=data$y-data$SB, p.bkg=p.bkg, sigma=data$sigma, lambda=data$lambda)

  max.log <- pmax(psi.likelihood.bkg, psi.likelihood.signal)  # parallel max --> vector of max values
  #  avoiding inf values:
  psi.likelihood <- max.log + log(exp(psi.likelihood.bkg - max.log) + exp(psi.likelihood.signal - max.log))
  psi.likelihood <- -sum(psi.likelihood)

  return(psi.likelihood)
}



####################################################################################
# logProbabilityBkgR(y, knots.y, alpha, Phi, Gr)
#
# returns:      the negative log of G(r)-part of the likelihood.
# arguments
#   y:          signal values.
#   knots.y:    the spline knot y-values ('c' in Fischer et al.).
#   alpha:      scale parameter.
#   Phi:        basis matrix.
#   Gr:         low-r G(r) information.
logProbabilityBkgR <- function(y, knots.y, alpha, Phi, Gr){
	if(is.na(Gr$type1))
      psi.gr.r1 <- 0
    else if(Gr$type1=="gaussianNoise")
	  psi.gr.r1 <- logLikelihoodGrGauss(y=y, knots.y=knots.y, alpha=alpha, Phi=Phi, bkg.r=Gr$bkg.r,
	                                 sigma.r=Gr$sigma.r, matrix.FT=Gr$matrix.FT1)$likelihood						 
	else if(Gr$type1=="correlatedNoise")
	  psi.gr.r1 <- logLikelihoodGrCorr(knots.y=knots.y, Phi=Phi, bkg.r=Gr$bkg.r, 
	                                KG.inv=Gr$KG.inv, matrix.FT=Gr$matrix.FT1)$likelihood
    if(is.na(Gr$type2))
      psi.gr.r2 <- 0
	else if(Gr$type2=="secondDeriv")
	  psi.gr.r2 <- logPriorBkgRSmooth(bkg.r=Gr$matrix.FT2 %*% (Phi %*% t(t(knots.y))), D=Gr$D)$likelihood
	else if(Gr$type2=="gaussianProcess")
      psi.gr.r2 <- logPriorBkgRGP(bkg.r=Gr$matrix.FT2 %*% (Phi %*% t(t(knots.y))), covMatrix=Gr$covMatrix)$likelihood
  return(psi.gr.r1 + psi.gr.r2)
}


########################################
# logPriorF(knots.x, knots.y)
#
# returns:     the negative log of the prior (due to Fischer et al.)
# arguments:    
#   knots.x:   the spline knot x-values ('xi' in Fischer et al.)
#   knots.y:   the spline knot y-values ('c' in Fischer et al.)
logPriorF <- function(knots.x, knots.y, Hessian=FALSE){
  E <- length(knots.y)
  if (length(knots.x) != E) {
    stop("Inconsistent lengths of parameters knots.y and knots.x!")
  }

  D <- DMatrix(knots.x=knots.x)
  detD <- D$det
  D <- D$matrix
  cDc <- as.vector(t(knots.y) %*% D %*% t(t(knots.y)))
  # Prior: separate terms according to which variables contribute
  prior.E  <- (0.5 * E * log(pi) - log(gamma(0.5 * E))) - sum(log(1:E)) 
  # "Depends" on knots.x, but constant as long as knots.x OK.
  prior.E.xi    <- -0.5 * log(detD)
  prior.E.xi.cc <- 0.5 * E * log(cDc)
  prior <- prior.E + prior.E.xi + prior.E.xi.cc
 
  hess <- NA
  if(Hessian==TRUE)
    hess <- E/2 * (2*(D / cDc) - (4 * D %*% t(t(knots.y)) %*% t(knots.y) %*% D) / (cDc ^ 2) )
  
  return (list(prior=prior, grad=(t(knots.y) %*% D) * E / cDc, hess=hess))
  
}

#####################################################################################
# logPriorBkgRSmooth(bkg.r, D, Hessian, Phi, matrix.FT, knots.y)
#
# returns:       prior in the r-space
# arguments
logPriorBkgRSmooth <- function(bkg.r, D, Hessian=FALSE, Phi=NA, matrix.FT=NA, knots.y=NA){
  
  E <- length(bkg.r)	
  cDc <- as.vector(t(bkg.r) %*% D %*% t(t(bkg.r)))
  # Prior: separate terms according to which variables contribute
  prior.E  <- 0.5 * E * log(pi) - log(gamma(0.5 * E)) - sum(log(1:E))
  # "Depends" on knots.x, but constant as long as knots.x OK.
#  prior.E.xi    <- -0.5 * log(detD)
  prior.E.xi    <- 0
  prior.E.xi.cc <- 0.5 * E * log(cDc)
  prior <- prior.E + prior.E.xi + prior.E.xi.cc
  
  hess <- NA
  if(Hessian==TRUE){
    D <- t(Phi) %*% t(matrix.FT) %*% D %*% matrix.FT %*% Phi
    hess <- E/2 * (2*(D / cDc) - (4 * D %*% t(t(knots.y)) %*% t(knots.y) %*% D) / (cDc ^ 2) )
  }
  return (list(likelihood=prior, hess=hess))
}

#####################################################################################
# logPriorBkgRGP(bkg.r, covMatrix, Hessian)
#
# returns:       prior in the r-space due to covariance matrix
# arguments
logPriorBkgRGP <- function(bkg.r, covMatrix, Hessian=FALSE){
  
  N <- length(bkg.r)
  f <- covMatrix$factor
  term.volume <- 0.5*N*log(2*pi)
  term.det <- 0.5*log(covMatrix$det) - 0.5*N*log(f)
  term.bkg <- f*0.5*t(bkg.r) %*% covMatrix$inv %*% t(t(bkg.r))
  
  prior <-  term.volume + term.det + term.bkg
  
  hess <- NA
  if(Hessian==TRUE)
	hess <- f*covMatrix$inv
 
  return (list(likelihood=prior, hess=hess))
  
}


#####################################################################################
# logLikelihoodBkg(knots.y, y, Phi, p.bkg, sigma) 
#
# returns:       the background-only contribution to the likelihood
# arguments
#   knots.y:     'c' in Fischer et al.: vector of spline knot y-values
#   y:           vector of datapoints
#   Phi:         spline matrix taking 'c' into 'y'
#   p.bkg:       the probability that a single pixel contains "only" background.
#   sigma:       vector of experimental uncertainties
logLikelihoodBkg <- function(knots.y, y, Phi, p.bkg, sigma) {
  deviation <- (y - Phi %*% t(t(knots.y)))  # 'y-bkg'_i should be eps_i
  deviation.norm <- deviation / sigma
  likelihood <- log(p.bkg) - 0.5 * log(2 * pi) - log(sigma) - 0.5 * deviation.norm ^ 2
  likelihood
}


logLikelihoodBkgAnalyt <- function(pars, x, y, p.bkg, sigma) {
  deviation <- (y - bkg.analyt(pars,x))  # 'y-bkg'_i should be eps_i
  deviation.norm <- deviation / sigma
  likelihood <- log(p.bkg) - 0.5 * log(2 * pi) - log(sigma) - 0.5 * deviation.norm ^ 2
  likelihood
}
#####################################################################################
# logLikelihoodSignal(knots.y, y, Phi, p.bkg, sigma, lambda)
#
# returns:       the signal-containing contribution to the likelihood
# arguments
#   knots.y:     'c' in Fischer et al.: vector of spline knot y-values
#   y:           vector of datapoints; 'd' at Fisher et al.
#   Phi:         spline matrix taking 'c' into 'y'
#   p.bkg:       The probability that a single pixel contains "only" background.
#   sigma:       vector of experimental uncertainties
#   lambda:      mean signal magnitude for signal-containing pixels. Either
#                vector of length y or a scalar value
logLikelihoodSignal <- function(knots.y, y, Phi, p.bkg, sigma, lambda) {
  rho <- sigma / lambda
  z <- (y - (Phi %*% t(t(knots.y)))) / lambda
  qq <- z / rho - rho
  likelihood <- log(1 - p.bkg) - log(lambda) + pnorm(log.p=TRUE, q=qq) - z + 0.5 * (rho ^ 2)
  likelihood
}

logLikelihoodSignalAnalyt <- function(pars, x, y, p.bkg, sigma, lambda) {
  rho <- sigma / lambda
  z <- (y - bkg.analyt(pars,x)) / lambda
  qq <- z / rho - rho
  likelihood <- log(1 - p.bkg) - log(lambda) + pnorm(log.p=TRUE, q=qq) - z + 0.5 * (rho ^ 2)
  likelihood
}
#####################################################################################
# logLikelihoodGrCorr(knots.y, Phi, bkg.r, KG.inv, matrix.FT, Hessian)
#
# returns:      the 'G(r)=-4*Pi*rho*r restriction' contribution to the likelihood.
# arguments
logLikelihoodGrCorr <- function(knots.y, Phi, bkg.r, KG.inv, matrix.FT, Hessian=FALSE){  

  M.Phi <- matrix.FT %*% Phi
  KG.inv.M.Phi <- KG.inv %*% M.Phi
  mu <- t(KG.inv.M.Phi) %*% bkg.r  # KG.inv is a symmetric matrix
  J <- t(KG.inv.M.Phi) %*% M.Phi
  psi.gr <- 0.5 * (t(knots.y) %*% J) %*% t(t(knots.y)) - t(mu) %*% t(t(knots.y))
  psi.gr <- sum(psi.gr)

  hess <- NA
  if(Hessian==TRUE)
	hess <- J
 
  return (list(likelihood=psi.gr, hess=hess))
}

#####################################################################################
#logLikelihoodGrGauss(y, knots.y, alpha, Phi, bkg.r, sigma.r, matrix.FT, Hessian=FALSE)
#
# returns:      the 'G(r)=-4*Pi*rho*r restriction' contribution to the likelihood.
# arguments  
logLikelihoodGrGauss <- function(y, knots.y, alpha, Phi, bkg.r, sigma.r, matrix.FT, Hessian=FALSE){  
  ############
  # y = SIGNAL - SB
  deviation  <- bkg.r/alpha + matrix.FT %*% ( (1-1/alpha)*y - Phi %*% t(t(knots.y)) )
  
  deviation.norm <- deviation / sigma.r
  psi.gr <- 0.5 * log(2 * pi) + log(sigma.r) + 0.5 * deviation.norm^2
  psi.gr <- sum(psi.gr)
   
  hess <- NA
  if(Hessian==TRUE){
    Mprime <- matrix.FT%*%Phi / (sqrt(2)*sigma.r)
	hess <- 2*t(Mprime) %*% Mprime
  }
  return (list(likelihood=psi.gr, hess=hess))
} 
 



#########################################################
# basisSpline(x, knots.x, knots.i, deriv)
#
# returns:       the i'th spline basis function (i.e. that is nonzero (=1) 
#                only in the ith knot) or its derivative at points x.
# arguments     
#    x:          x-values where we evaluate the basis functions.
#    knots.x:    knot positions
#    knots.i:    which basis function to return
#    deriv:      which derivative of the spline function (0 to 3)
basisSpline <- function(x, knots.x, knots.i, deriv=0){
  y <- 0*knots.x
  y[knots.i] <- 1
  basisSpline <- splinefun(x=knots.x, y=y, method="natural")
  return(basisSpline(x=x, deriv=deriv))
}
 
#########################################################
# basisMatrix(x, knots.x,deriv)
#
# returns:       matrix whose columns are spline basis functions 
#                values at points x; if an arbitrary function is
#                known only at points knots.x, multiplying it by
#                basisMatrix results in function spline approximations 
#                at points x
# arguments
#   x:           x-values where we evaluate the basis functions.
#   knots.x:     knot positions
#   deriv:       which derivative of the spline function?
basisMatrix <- function(x, knots.x, deriv=0) {
  E <- length(knots.x)
  N <- length(x)
  bM <- matrix(nrow=N, ncol=E)
  for (i in 1:E) {
    bM[, i] <- basisSpline(x=x, knots.x=knots.x, knots.i=i, deriv=deriv)
  }
  return(bM)
}
 

########################################################################################
# DMatrix(knots.x, only.trace, robust.factor)
#
# returns:          the second derivative overlap matrix ('D' in Fischer et al.
#                   (1999)) of the spline basis functions with knots at 'knots.x'.
# arguments    
#   knots.x:        spline knot positions
#   onlyTrace:      if true, returns trace(D) (a single number) instead of D.
#   robustFactor:   We add a constant to all eigenvalues to promote stability;
#                   this is the ratio of that constant to the smallest
#                   eigenvalue which is SUPPOSED to be nonzero.
DMatrix <- function(knots.x, onlyTrace=FALSE, robustFactor=1e-12) {
 # cat("Calculating basis matrix...")
  N <- length(knots.x)
  ddM <- basisMatrix(x=knots.x, knots.x=knots.x, deriv=2)
 # cat(" done!\n")
  x.mat <- matrix(rep(knots.x, N), nrow=N)
  dy <- apply(ddM, 2, diff)
  dx <- apply(x.mat, 2, diff)
  # slope[i, j] is the slope of the i'th segment of the j'th basis function;
  # intercept[i, j] is its y-intercept
  slope <- dy / dx
  intercept <- ddM[-1, ] - slope * x.mat[-1, ]
  d.x3 <- diff(knots.x ^ 3)
  d.x2 <- diff(knots.x ^ 2)
  # The for-loop code below reduces to this if we only care about diagonal
  # elements (which is the case for computing the trace).
  if (onlyTrace) {
    return (sum((slope ^ 2) * d.x3 / 3.0 + (slope * intercept) * d.x2 + 
	    intercept ^ 2 * diff(knots.x)))
  }
  # If we're this far, we need to compute the whole matrix.
  DD <- matrix(0, nrow=N, ncol=N)
  
  # computing basis function support:
  N2 <- floor(N/2)
  phi.N2 <- basisSpline(x=knots.x, knots.x=knots.x, knots.i=N2, deriv=2)
  max.phi <- max(phi.N2)
  x3.factor <- (max(knots.x)^3-min(knots.x)^3)/3
  supp <- min(length(which(abs(phi.N2)>1e-9*max.phi/x3.factor)), N-1)  # only i-supp to i+supp matters
  
  # cat("Calculating overlap intergals...\n")
  for (i in 1:N) {
    #    if (i %% 100 == 0)
    #	  cat("...x=", knots.x[i], "\n")
	relevant <- min(i+supp, N)
	phi.i.supp.min <- max(1, i-ceiling(supp/2))
	phi.i.supp.max <- min(N-1, i+floor(supp/2)) 
    for (j in i:relevant) {
      phi.j.supp.min <- max(1, j-ceiling(supp/2))
	  phi.j.supp.max <- min(N-1, j+floor(supp/2)) 
	  phi.ij.supp <- max(phi.i.supp.min, phi.j.supp.min):min(phi.i.supp.max, phi.j.supp.max)
      prod.ss <- slope[phi.ij.supp, i] * slope[phi.ij.supp, j] / 3.0
      prod.si <- 0.5 * (slope[phi.ij.supp, i] * intercept[phi.ij.supp, j] + slope[phi.ij.supp, j] * intercept[phi.ij.supp, i])
      prod.ii <- intercept[phi.ij.supp, i] * intercept[phi.ij.supp, j]
      integrals <- prod.ss * d.x3[phi.ij.supp] + prod.si * d.x2[phi.ij.supp] + prod.ii * diff(knots.x)[phi.ij.supp]
      DD[i, j] <- DD[j, i] <- sum(integrals)
    }
  }
  #  cat("...done!\n")
  #  cat("Calculating D-matrix eigenvalues...")
  # "modified" determinant (i.e., product of highest (N-2) eigenvalues)
  eigenvals <- eigen(x=DD, symmetric=TRUE, only.values=TRUE)$values
  # cat(" done!\n")
  detD <- prod(rev(eigenvals)[-(1:2)])

  # Make sure we avoid negative eigenvalues!
  lowest.true.eigenval <- sort(eigenvals)[3]

  return (list(matrix=DD + lowest.true.eigenval * robustFactor * diag(nrow(DD)), det=detD))
}

#########################################################
# get.bkg(x, knots.x, knots.y)
#
# returns:    natural spline approximation of function knots.y
# arguments
#   x:        points at which function should be approximated 
#   knots.x:  knot position
#   knots.y:  knot values
get.bkg <- function(x, knots.x, knots.y){
  bM <- basisMatrix(x=x, knots.x=knots.x)
  get.bkg <- c(bM %*% t(t(knots.y)))
  return(get.bkg)
}


bkg.analyt <- function(pars, x){
  bkg.analyt <- pars[1]*exp(-pars[2]*x)*x^pars[3] + pars[4]/((x-pars[5])^2+pars[6]^2)
}
