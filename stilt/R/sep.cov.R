sep.cov <-
function(Theta.mat, t.vec, rho, kappa, phi.vec, zeta) {

  
# PRELIMINARIES #!+
n.par <- length(t.vec)
p.par <- dim(Theta.mat)[1]
m.par <- dim(Theta.mat)[2]
if (rho >= 1) stop('***ERROR** rho should be less than 1!')
if (rho < 0) stop('***ERROR** rho should be non-negative!')
if (zeta < 0) stop('***ERROR*** nugget should be non-negative')
if (any(phi.vec <= 0)) stop ('***ERROR*** All phi pars should be positive')
if (kappa < 0) stop('***ERROR*** Kappa should be non-negative')


# TIME COVARIANCE #!+
diff.mat    <- abs(outer(t.vec, t.vec, FUN='-'))
Sigma.t.mat <- (rho^diff.mat)/(1-rho^2) 

  
# PARAMETER COVARIANCE  #!+
# Range matrix corresponding to Theta.mat 
Phi.mat         <- matrix(phi.vec, nrow=p.par, ncol=m.par, byrow=TRUE)
# Parameter matrix scaled by phi 
Theta.mat.sc    <- Theta.mat/Phi.mat 
dist.mat        <- rdist(Theta.mat.sc) #Sigma term 
Sigma.theta.mat <- kappa*exp(-dist.mat^2) 
Sigma.theta.mat <- Sigma.theta.mat + zeta*diag(nrow=p.par) 


# FORMAT OUTPUT #!+
Sigma.mats <- list(Sigma.t.mat=Sigma.t.mat, Sigma.theta.mat=Sigma.theta.mat)
Sigma.mats
}
