emul.lik <-
function(parvec, Y.mat, X.mat, t.vec, Theta.mat, n.par, p.par, fix.betas,
                     limits.lower=NULL, limits.upper=NULL, beta.vec=NULL) {

  
# PRELIMINARIES #!+
if ((fix.betas) && (is.null(beta.vec))) {
  stop("***ERROR*** Betas are fixed, yet the beta vector was not provided!\n")
}
if ( (!fix.betas) && (!is.null(beta.vec))) {
  cat("WARNING: 'beta.vec' argument is ignored\n")
}
  
# EXTRACT PARAMETERS #!+
rho       <- parvec[1] 
kappa     <- parvec[2] 
zeta      <- parvec[3] 
if (!fix.betas) { 
   beta.ind  <- names(parvec) == "beta" 
   beta.vec  <- parvec[beta.ind]
}
phi.ind   <- names(parvec) == "phi"
phi.vec   <- parvec[phi.ind] 

# SOME CHECKS #!+
if (kappa == 0 && zeta == 0) stop("***ERROR*** Kappa and zeta can't be both 0!\n")

# SET-UP MATRICES #!+
Sigma.mats          <- sep.cov(Theta.mat, t.vec, rho, kappa, phi.vec, zeta)

#browser()

Sigma.theta.inv.mat <- solve(Sigma.mats$Sigma.theta.mat) 
Sigma.t.inv.mat     <- solve(Sigma.mats$Sigma.t.mat) 
mu.vec              <- X.mat%*%as.matrix(beta.vec) 
vec.C               <- Y.mat - mu.vec 
C.mat               <- matrix(as.vector(vec.C), nrow=p.par, ncol=n.par)


# CALCULATE LIKELIHOOD #!+
# The new determinant function behaves fine for matrices with very large and close to 0
# numbers
T1.mat              <- Sigma.theta.inv.mat%*%C.mat%*%Sigma.t.inv.mat 
T2.mat              <- C.mat*T1.mat 
Term1               <- -0.5*sum(T2.mat)
# Det10, Det20: 1st element is ln(modulus of the determinant), the 2nd is the sign of the determinant
Det10               <- unname(unlist(determinant(Sigma.mats$Sigma.t.mat, logarithm=TRUE)))
Det1                <- p.par*Det10[1] 
Det20               <- unname(unlist(determinant(Sigma.mats$Sigma.theta.mat,
                                                 logarithm=TRUE))) 
Det2                <- n.par*Det20[1] 
Term2               <- -0.5*(Det1+Det2) 
Term3               <- -0.5*n.par*p.par*log(2*pi) 
llik                 <- Term1 + Term2 + Term3 


# CHECK THAT DETERMINANTS ARE POSITIVE #!+
if ((Det10[2] < 0) || (Det20[2] < 0)) {
  stop("***ERROR*** Covariance matrix determinant(s) is/are negative!")
}


# IF OUTSIDE LIMITS ASSUME NEGATIVE INFINITI #!+
#!+
if (!is.null(limits.lower)) { 
  if (any(parvec < limits.lower)) {
    llik <- -Inf #!+
  }
}
#!+
if (!is.null(limits.upper)) { 
  if (any(parvec > limits.upper)) {
    llik <- -Inf #!+
  }
}

#cat('Log-likelihood=', llik, '\n')

#!+
llik 
}
