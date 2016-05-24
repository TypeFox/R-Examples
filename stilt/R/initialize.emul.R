initialize.emul <-
function(mpars, moutput, par.reg, time.reg, kappa0, zeta0) {

  
# PRELIMINARIES #!+
  m.par     <- dim(mpars$par)[1] #Number of parameters 
  p.par     <- dim(mpars$par)[2] #Number of ensemble members 
  n.par     <- dim(moutput$out)[1] #Length of time dimension 
  nreg      <- 1 + sum(par.reg) + sum(time.reg) #Number of regressors (including the constant term) 

  
# OBTAIN DATA VECTOR, PARAMETER MATRIX, DESIGN MATRIX AND COVARIATES MATRIX #!+
  dmat      <- design.mat(mpars, moutput, par.reg, time.reg) 
  Theta.mat <- dmat$Theta.mat 
  t.vec     <- moutput$t

  
# INITIAL GUESS FOR BETA PARAMETERS #!+
# Estimated beta values are in beta.est$coefficients
  X.mat.df  <- as.data.frame(dmat$X.mat) 
  beta.est  <- lm(dmat$Y.mat ~ 0 + ., data=X.mat.df, offset=NULL) 
  beta.vec  <- unname(beta.est$coefficients)

# INITIAL GUESS FOR RHO #!+
  rho       <- 0.9 

# INITIAL GUESS FOR RANGE PARAMETERS #!+
# (Half of the parameter ranges covered by the ensemble)
  par.min   <- apply(Theta.mat, 2, min) #Min of each parameter 
  par.max   <- apply(Theta.mat, 2, max) #Max of each parameter 
  phi.vec   <- (par.max-par.min)/2 



# CALCULATE MEAN VECTOR AND VECC #!+
  mu.vec    <- dmat$X.mat%*%as.matrix(beta.vec)
  vecC      <- dmat$Y.mat - mu.vec #w


# CALCULATE COVARIANCE MATRIX #!+
  Sigma.mats <- sep.cov(Theta.mat, t.vec, rho, kappa0, phi.vec, zeta0)
  

# CONSTRUCT THE EMULATOR #!+
  init.emul <- list(Theta.mat=Theta.mat, t.vec=t.vec, Y.mat=dmat$Y.mat, X.mat=dmat$X.mat,
                    beta.vec=beta.vec, kappa=kappa0, phi.vec=phi.vec, zeta=zeta0, n=n.par,
                    rho=rho, p=p.par, vecC=vecC, par.reg=par.reg, time.reg=time.reg)
  init.emul
}
