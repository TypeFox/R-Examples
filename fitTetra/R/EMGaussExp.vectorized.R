EMGaussExp.vectorized <-
function(y,psi) {
  ng <- length(psi$mu)
  ye <- rep(y,each=ng)
  dens.norm <- matrix(dnorm(ye, psi$mu, psi$sigma), ncol=ng, byrow=T)
  Zp <- dens.norm %*% diag(psi$p)
  Zs <- dens.norm %*% psi$p
  Z <- t(scale(t(Zp), center=F, scale=Zs))
  attributes(Z)$"scaled:scale" <- NULL
  Z    # Z is matrix of posterior probabilities with as many row as length(y) and ng columns
}
