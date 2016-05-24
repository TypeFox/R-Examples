FLLat.BIC <- function(Y,J=min(15,floor(ncol(Y)/2)),B="pc",thresh=10^(-4),
                      maxiter=100,maxiter.B=1,maxiter.T=1) {

  ## Error checking parameters.
  CheckPars(Y=Y,J=J,B=B,thresh=thresh,maxiter=maxiter,maxiter.B=maxiter.B,
            maxiter.T=maxiter.T)

  ## Initializing Beta.
  if (is.matrix(B)) {
    in.B <- B
  } else if (B=="pc") {
    Y.cen <- scale(Y,scale=FALSE)
    Y.v <- svd(Y.cen)$v
    in.B <- Y.cen%*%(Y.v[,1:J,drop=FALSE])
  } else if (B=="rand") {
    in.B <- Y[,sample(ncol(Y),J),drop=F]
  }

  ## Setting up parameters.
  n <- length(Y)
  n.lam0s <- 5
  n.alphas <- 5
  alphas <- seq(0.1,0.9,len=n.alphas)
  in.lam0 <- max(abs(in.B))*seq(3,1,len=n.alphas)

  ## Results.
  lam0.max <- rep(0,n.alphas)
  bics <- lam0s <- matrix(0,nrow=n.alphas,ncol=n.lam0s)

  ## Generating grid of lams (alphas in rows, lam0s in columns).
  for (i in 1:n.alphas) {
    lam0.max[i] <- Max.Lam0(Y,J,in.B,in.lam0[i],alphas[i],thresh,maxiter,maxiter.B,
                            maxiter.T)
    lam0s[i,] <- seq(0,lam0.max[i],len=n.lam0s+2)[-c(1,n.lam0s+2)]
  }
  
  ## Calculating BICs.
  for (j in 1:n.lam0s) {
    for (i in 1:n.alphas) {
      bic.est <- FLLat(Y,J,in.B,lam0s[i,j]*alphas[i],
                       lam0s[i,j]*(1-alphas[i]),thresh,maxiter,maxiter.B,
                       maxiter.T)
      bics[i,j] <- bic.est$bic
      if (all(i==1,j==1)) {
        opt.est <- bic.est
        cur.bic <- bics[i,j]
        opt.i <- i; opt.j <- j
      } else if (bics[i,j]<=cur.bic) {
        opt.est <- bic.est
        cur.bic <- bics[i,j]
        opt.i <- i; opt.j <- j
      }
    }
  }

  ## Optimum lam0 and alpha.
  opt.lam0 <- lam0s[opt.i,opt.j]; opt.alpha <- alphas[opt.i]
  
  return(list("lam0"=opt.lam0,"alpha"=opt.alpha,"lam1"=opt.lam0*opt.alpha,
              "lam2"=opt.lam0*(1-opt.alpha),"opt.FLLat"=opt.est))
  
}
