pfa <-
function(X, option="data", m=2, initial.communality ="R2", crit=0.001, verbose = FALSE)
{
S <- switch(option, cor = X, cov = X, data = cor(X), stop("pfa: invalid value for parameter option")) 
p <- ncol(X)
R <- diag(1/sqrt(diag(S)))%*%S%*%diag(1/sqrt(diag(S)))
if(initial.communality=="R2")
  #comu <- diag(1-1/solve(R))*diag(S)
  comu <- 1-1/diag(solve(R))
if(initial.communality=="maxcor") {
  Rn <- abs(R)
  diag(Rn) <- 0
  comu <- apply(Rn,1,max)
}
cat("Initial communalities\n")
print(comu)
L <- diag(svd(S)$d)
V <- svd(S)$u
A <- as.matrix(V[,1:m])%*%as.matrix(sqrt(L[1:m,1:m]))
eps <- 1
epsh <- c(eps)
niter <- 1
C <- S 
while (eps > crit) {
   diag(C) <- comu 
   L <-diag(svd(C)$d)
   V <- svd(C)$u
   A <- as.matrix(V[,1:m])%*%as.matrix(sqrt(L[1:m,1:m]))
   comu.r <- diag(A%*%t(A))
   corr <-(1:length(comu.r)) [comu.r > diag(S)]
   comu.r[corr] <- diag(S)[corr]
   eps <- max(abs(comu -  comu.r))
   epsh <- c(epsh, eps)
   comu <- comu.r
   niter <- niter + 1
}
if(verbose) {
  ithist <- cbind(1:niter,epsh)
  cat("iteration history\n")
  print(ithist)
}
La <- A
Psi <- diag(diag(S)- comu)
cat("Final communalities\n")
print(comu)
cat(niter," iterations till convergence\n")
cat("Specific variances:\n")
print(diag(Psi))
Shat <- (La%*%t(La) + Psi)
cat("Variance explained by each factor\n")
Dv <- diag(t(La)%*%La)
print(Dv)
cat("Loadings:\n")
print(La)
Res <- S - La%*%t(La) + Psi
if(option=="data") {
  Xs <- scale(X)
  Fs <- Xs%*%solve(R)%*%La
} else Fs <- NA
list(Res = Res, Psi = Psi, La = La, Shat = Shat, Fs = Fs)
}

