
lisrel <- function(model,p,X=NULL,muX=NULL,varX=NULL,...) {
  pp <- modelPar(model,p)
  mom <- moments(model,p)
  A <- t(index(model)$M)
  J <- index(model)$J ## Observed var. selection matrix

  eta.idx <- match(latent(model),vars(model))
  obs.idx <- match(manifest(model),vars(model))
  exo.idx <- match(exogenous(model),vars(model))
  y <- setdiff(manifest(model), exogenous(model))
  y.idx <- match(y, vars(model))

##  Jy <- Jx <- Jeta <- I <- diag(length(vars(model)))
##  if (length(eta.idx)>0)
##    J[eta.idx,eta.idx] <- 0; J <- J[-eta.idx,]
##  Jeta[obs.idx,obs.idx] <- 0; Jeta <- J[-obs.idx,]

  A <- t(mom$A)
  Lambda <- A[y.idx,eta.idx,drop=FALSE]
  K <- A[y.idx,exo.idx,drop=FALSE]
  B <- A[eta.idx,eta.idx,drop=FALSE]
  I <- diag(nrow=nrow(B))
  Gamma <- A[eta.idx,exo.idx,drop=FALSE]
  V <- mom$P
  Psi <- V[eta.idx,eta.idx] ## Residual variance

  Theta <- V[y.idx,y.idx] ## -
  IBi <- if (ncol(I)>0) solve(I-B) else I
  LIBi <- Lambda%*%IBi
  Phi <- LIBi%*%Gamma + K

  Veta.x <- IBi%*%Psi%*%IBi  ## Variance of eta given x
  COVetay.x <- Veta.x%*%t(Lambda) ## Covariance of eta,y given x
##  Vy.x <- Lambda%*%COVetay.x + Theta ## Omega
  Vy.x <- LIBi%*%Psi%*%t(LIBi) + Theta

  if (!is.null(X)) {
    Ey.x <- t(apply(as.matrix(X)%*% t(LIBi%*%Gamma + K),1,function(x) x + mom$v[y.idx]))
  } else Ey.x <- NULL

  Sigma <- mom$Cfull
  CV <- COVetay.x%*%Vy.x

##  Sigma <- Vy.x + Phi%*%varX%*%t(Phi)

  return(list(Lambda=Lambda, K=K, B=B, I=I, Gamma=Gamma, Psi=Psi, Theta=Theta, IBi=IBi, LIBi=LIBi, Phi=Phi,
         Vy.x=Vy.x, Veta.x=Veta.x, COVetay.x=COVetay.x, CV=CV, Ey.x=Ey.x))
}
