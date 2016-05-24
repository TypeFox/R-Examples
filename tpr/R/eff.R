## #### one-step update from tpr using sparse matrix operations
## #### get H and G for 1step update
## eff1step <- function(y, delta, family, x, z, tis, fit, var=TRUE, phi=0.95) {
##   ## fit: a tpr fit
##   n <- length(y)
##   nt <- length(tis)
##   p <- ncol(x)
##   q <- ncol(z)

##   ##phi <- 0.95 ## fixed for now
  
##   P <- nt * p + q ## the dimension of BETA
##   BETA <- c(as.vector(t(fit$alpha)), fit$theta)
##   H <- Matrix(0, P, P)
##   G <- Matrix(0, P, 1)
##   ##  GiGiT <- Matrix(0, P, P)
##   for (i in 1:n) {
##     ##cat("i = ", i, "\n")
##     DLTi <- interpprev(delta[[i]], tis)
##     if (sum(DLTi == 0)) next
##     Xi <- getXi(x[i,], z[i,], nt)
##     Yi <- interpprev(y[[i]], tis)     
##     HiGi <- getHiGi(Yi, Xi, DLTi, family, BETA, phi)
##     H <- H + HiGi$Hi
##     G <- G + HiGi$Gi
##     ## if (var) GiGiT <- GiGiT + tcrossprod(HiGi$Gi)
##   }
##   chg <- solve(H, G)
##   BETA <- BETA + chg
  
##   alpha <- matrix(as.vector(BETA[1:(nt * p)]), nrow=nt, ncol=p, byrow=TRUE)
##   theta <- BETA[-(1:(nt * p))]
##   valpha <- vtheta <- NULL
##   ## compute variance
##   if (var) {
##     H <- Matrix(0, P, P)
##     G <- Matrix(0, P, 1)
##     GiGiT <- Matrix(0, P, P)
##     for (i in 1:n) {
##       DLTi <- interpprev(delta[[i]], tis)
##       if (sum(DLTi == 0)) next
##       Xi <- getXi(x[i,], z[i,], nt)
##       Yi <- interpprev(y[[i]], tis)     
##       HiGi <- getHiGi(Yi, Xi, DLTi, family, BETA, phi)
##       H <- H + HiGi$Hi
##       GiGiT <- GiGiT + tcrossprod(HiGi$Gi)
##     }
##     vBETA <- diag(solve(H, GiGiT) %*% solve(H))
##     valpha <- matrix(vBETA[1:(nt * p)], nrow=nt, ncol=p, byrow=TRUE)
##     vtheta <- vBETA[-(1:(nt * p))]
##   }
##   list(tis=tis, alpha=alpha, theta=theta, valpha=valpha, vtheta=vtheta)
## }

## #### This function gets Xi, watch for the dimension:
## #### ETA = crossprod(Xi, BETA), so Xi is P by NT (BETA is P by 1)
## getXi <- function(xi, zi, NTi) {
##   ## xi: p by 1 matrix, Xmat[i,]
##   ## zi: q by 1 matrix, Zmat[i,]
##   ## NTi: the number of time points for subject i
##   Xi <- kronecker(Diagonal(NTi), as.matrix(xi)) ## this is much faster
##   ##Xi <- bdiag(rep(list(as.matrix(xi)), NTi))
##   if (length(zi) > 0) Xi <- rBind(Xi, rep(zi, NTi))
##   Xi
## }

## getQi <- function(DLTi, phi) {
##   nt <- length(DLTi)
##   Qi <- Matrix(0, nt, nt)
    
##   ## the first nonzero delta, assuming consecutibity
##   idx <- (1:nt)[DLTi == 1][1] 
  
##   NTi <- sum(DLTi)
##   if (NTi == 0) return(Qi)

##   if (NTi == 1) {
##     Qi[idx, idx] <- 1
##   }
##   else {
##     Qi[idx -1 + 1:NTi, idx - 1 + 1:NTi] <- Diagonal(x=c(1, rep(1 + phi*phi, NTi - 2), 1))
##     Qi[idx - 1 + 2:NTi, idx - 1 + 1:(NTi - 1)] <-
##       Qi[idx - 1 + 2:NTi, idx - 1 + 1:(NTi - 1)] + Diagonal(x = rep(-phi, NTi - 1))
##     Qi[idx - 1 + 1:(NTi - 1), idx - 1 + 2:NTi] <-
##       Qi[idx - 1 + 1:(NTi - 1), idx - 1 + 2:NTi] + Diagonal(x = rep(-phi, NTi - 1))
##   }
##   Qi
## }

## #### the cholesky decomposition of Qi 
## getLQi <- function(NTi, phi) {
##   if (NTi == 1) return(Diagonal(1))
##   LQi <- Diagonal(x=c(rep(1, NTi - 1), sqrt(1 - phi*phi)))
##   ##ind <- outer(1:NTi, 1:NTi, function(x, y) x - y == 1)
##   ##LQi[ind] <- -phi
##   LQi[2:NTi, 1:(NTi - 1)] <- LQi[2:NTi, 1:(NTi - 1)] -
##     Diagonal(x = rep(phi, NTi - 1))
##   LQi
## }

## getHiGi <- function(Yi, Xi, DLTi, family, BETA, phi) {
##   ## good <- DLTi == 1
##   ETA <- as.vector(crossprod(Xi, BETA))
##   MU <- family$linkinv(ETA)
##   MU.ETA <- family$mu.eta(ETA)
##   Vroot <- sqrt(family$variance(MU))
##   Di <- Xi %*% Diagonal(x = MU.ETA * DLTi)
##   Qi <- getQi(DLTi, phi)
##   OMEGAi <- Diagonal(x = Vroot) %*% Qi %*% Diagonal(x = Vroot)
##   Ei <- Matrix(Yi - MU, ncol = 1)
##   DLT <- Diagonal(x = DLTi)
##   ## output
##   Gi = Di %*% crossprod(OMEGAi,  Ei)
##   Hi = Di %*% tcrossprod(OMEGAi, Di)
##   list(Hi = Hi, Gi = Gi)
## }
