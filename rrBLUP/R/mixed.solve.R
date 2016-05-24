mixed.solve <- function (y, Z = NULL, K = NULL, X = NULL, method = "REML", bounds = c(1e-09,1e+09), SE = FALSE, return.Hinv = FALSE) {
pi <- 3.14159
n <- length(y)
y <- matrix(y,n,1)

not.NA <- which(!is.na(y))

if (is.null(X)) {
p <- 1
X <- matrix(rep(1,n),n,1)
}
p <- ncol(X)
if (is.null(p)) {
p <- 1
X <- matrix(X,length(X),1)
}
if (is.null(Z)) {
Z <- diag(n)
}
m <- ncol(Z)
if (is.null(m)) {
m <- 1
Z <- matrix(Z,length(Z),1)
}
stopifnot(nrow(Z) == n)
stopifnot(nrow(X) == n)
if (!is.null(K)) {
	stopifnot(nrow(K) == m)
	stopifnot(ncol(K) == m)
}

Z <- as.matrix(Z[not.NA,])
X <- as.matrix(X[not.NA,])
n <- length(not.NA)
y <- matrix(y[not.NA],n,1)

XtX <- crossprod(X, X)
rank.X <- qr(XtX)$rank
if (rank.X < p) {stop("X not full rank")}
XtXinv <- solve(XtX)
S <- diag(n) - tcrossprod(X%*%XtXinv,X)
if (n <= m + p) {
  spectral.method <- "eigen"
} else {
  spectral.method <- "cholesky"
  if (!is.null(K)) {
    diag(K) <- diag(K) + 1e-6
  	B <- try(chol(K),silent=TRUE)
  	if (class(B)=="try-error") {stop("K not positive semi-definite.")}
  } # if is.null
} 

if (spectral.method=="cholesky") {
if (is.null(K)) {
	ZBt <- Z
} else {
	ZBt <- tcrossprod(Z,B) 
}
svd.ZBt <- svd(ZBt,nu=n)
U <- svd.ZBt$u
phi <- c(svd.ZBt$d^2,rep(0,n-m))
SZBt <- S %*% ZBt
svd.SZBt <- try(svd(SZBt),silent=TRUE)
if (class(svd.SZBt)=="try-error") {svd.SZBt <- svd(SZBt+matrix(1e-10,nrow=nrow(SZBt),ncol=ncol(SZBt)))}
QR <- qr(cbind(X,svd.SZBt$u))
Q <- qr.Q(QR,complete=TRUE)[,(p+1):n]
R <- qr.R(QR)[p+1:m,p+1:m]
ans <- try(solve(t(R^2), svd.SZBt$d^2),silent=TRUE)
if (class(ans)=="try-error") {
    spectral.method <- "eigen"
} else {
    theta <- c(ans,rep(0, n - p - m))
}
}

if (spectral.method=="eigen") {
offset <- sqrt(n)
if (is.null(K)) {
	Hb <- tcrossprod(Z,Z) + offset*diag(n)
} else {
	Hb <- tcrossprod(Z%*%K,Z) + offset*diag(n)
}
Hb.system <- eigen(Hb, symmetric = TRUE)
phi <- Hb.system$values - offset
if (min(phi) < -1e-6) {stop("K not positive semi-definite.")}
U <- Hb.system$vectors
SHbS <- S %*% Hb %*% S
SHbS.system <- eigen(SHbS, symmetric = TRUE)
theta <- SHbS.system$values[1:(n - p)] - offset
Q <- SHbS.system$vectors[, 1:(n - p)]
}  
    
omega <- crossprod(Q, y)
omega.sq <- omega^2
if (method == "ML") {
f.ML <- function(lambda, n, theta, omega.sq, phi) {
 n * log(sum(omega.sq/(theta + lambda))) + sum(log(phi + lambda))
}
soln <- optimize(f.ML, interval = bounds, n, theta, omega.sq, phi)
lambda.opt <- soln$minimum
df <- n
} else {
f.REML <- function(lambda, n.p, theta, omega.sq) {
 n.p * log(sum(omega.sq/(theta + lambda))) + sum(log(theta + lambda))
}
soln <- optimize(f.REML, interval = bounds, n - p, theta, omega.sq)
lambda.opt <- soln$minimum
df <- n - p
} #if method
Vu.opt <- sum(omega.sq/(theta + lambda.opt))/df
Ve.opt <- lambda.opt * Vu.opt
Hinv <- U %*% (t(U)/(phi+lambda.opt))
W <- crossprod(X,Hinv%*%X)
beta <- array(solve(W,crossprod(X,Hinv%*%y)))
rownames(beta) <- colnames(X)

if (is.null(K)) {
	KZt <- t(Z)
} else {
	KZt <- tcrossprod(K,Z)
}
KZt.Hinv <- KZt %*% Hinv
u <- array(KZt.Hinv %*% (y - X%*%beta))

if (is.null(K)) {
    rownames(u) <- colnames(Z)
} else {
    rownames(u) <- rownames(K)
}

LL = -0.5 * (soln$objective + df + df * log(2 * pi/df))
if (!SE) {
  if (return.Hinv) {
    return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, u = u, LL = LL, Hinv = Hinv))
  } else {
    return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, u = u, LL = LL))
  }
} else {
  Winv <- solve(W)
  beta.SE <- array(sqrt(Vu.opt*diag(Winv)))
  rownames(beta.SE) <- rownames(beta)
  WW <- tcrossprod(KZt.Hinv,KZt)
  WWW <- KZt.Hinv%*%X
  if (is.null(K)) {
	u.SE <- array(sqrt(Vu.opt * (rep(1,m) - diag(WW) + diag(tcrossprod(WWW%*%Winv,WWW)))))	
  } else {
	u.SE <- array(sqrt(Vu.opt * (diag(K) - diag(WW) + diag(tcrossprod(WWW%*%Winv,WWW)))))
  }
  rownames(u.SE) <- rownames(u)
  
  if (return.Hinv) {
    return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, beta.SE = beta.SE, u = u, u.SE = u.SE, LL = LL, Hinv = Hinv))
  } else {
    return(list(Vu = Vu.opt, Ve = Ve.opt, beta = beta, beta.SE = beta.SE, u = u, u.SE = u.SE, LL = LL))
  }
}
}
