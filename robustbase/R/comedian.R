### -*- mode: R ; delete-old-versions: never -*-

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## http://www.r-project.org/Licenses/

### From package 'Biobase' (has only rowMedians + rowQ) / 'matrixStats'
### MM: all type checking now in C
## ---  TODO: implement  hasNA=NA ==> do check maybe differently than = TRUE
##  --> ../src/rowMedians.c +  ../src/rowMedians_TYPE-template.h
colMedians <- function(x, na.rm=FALSE, hasNA=TRUE, keep.names=TRUE)
  .Call(R_rowMedians, x, na.rm, hasNA, FALSE, keep.names)

rowMedians <- function(x, na.rm=FALSE, hasNA=TRUE, keep.names=TRUE)
  .Call(R_rowMedians, x, na.rm, hasNA, TRUE, keep.names)




### Maria Anna di Palma, without consistency factor 15.11.2014
### Fixes by Valentin Todorov
### Martin Maechler: added mad() consistency factor, 27.11.2014
###     new name, class; more compatible to 'covMcd'

covComed <- function (X, n.iter = 2, reweight = FALSE,
                      tolSolve = control$ tolSolve,# had 1e-10 hardwired {now 1e-14 default}
                      trace = control$ trace,
                      wgtFUN = control$ wgtFUN,
                      control = rrcov.control())
{
    ## ATTENTION ##
    ## Med(abs(X))^2=Med(X*X) only if the number of rows is odd
    d <- dim(X <- as.matrix(X))
    n <- d[1]
    p <- d[2]

    if(is.character(wgtFUN))
	wgtFUN <- .wgtFUN.covComed[[wgtFUN]](p=p, n=n, control)
    if(!is.function(wgtFUN))
	stop("'wgtFUN' must be a function or a string specifying such a function")

    madX <- apply(X, 2, mad)
    I.mad <- 1/madX
    rho <- I.mad * COM(X) * rep(I.mad, each = p)
    ## better than
    ## D <- diag(1/madX)
    ## rho <- D %*% COM(X) %*% t(D)

    U <- svd(rho, p, nv = 0L)$u
    ## DD <- diag(madX)
    ## Q <- DD %*% U
    ## invQ <- solve(Q) ##  == t(U) %*% D -- since  U is orthogonal!
    t.inv.Q <- I.mad * U # = t(solve(Q)) = t(t(U) * D) == t(D) U = D U
    Z <- X %*% t.inv.Q ## much faster than for (i in 1:n) Z[i,] <- invQ %*% X[i,]
    out <- comedian(rho, Z, X)

    ## Mahalanobis distance
    for(it in seq_len(n.iter))# allow n.iter = 0
        out <- comedian(out$S., out$Z, X)

    mm <- colMedians(out$Z)
    mx <- drop(out$Q %*% mm)
    ## MM: These are "raw" distances compared to covMcd()
    mah <- mahalanobis(X, mx, out$S., tol = tolSolve)

    ## compute weights
    weights <- wgtFUN(mah)
    covW <- cov.wt(X, wt=weights)[c("cov", "center", "n.obs")]
    covW$weights <-
	if(reweight) { ## above 'mah' = 'raw.mah' .. ==> allow another reweighting as in covMcd()
	    covW$raw.weights <- weights
	    covW$mah <- mahalanobis(X, covW$center, covW$cov, tol = tolSolve)
	    wgtFUN(mah)
	} else # no re-weighting
	    weights
    structure(class = "comed",
	      c(list(Z = out$Z, raw.cov = out$S., raw.center = mx, raw.mah = mah,
		     wgtFUN=wgtFUN),
		covW))
}


##' Martin Maechler's simple proposal for an *adaptive* cutoff
##' i.e., one which does *not* reject outliers in good samples asymptotically
.COM.adaptWgt.c <- function(n,p, eps = 0.2 / n^0.3) {
    ## default eps ==>  1-eps(n=100) ~= 0.95; 1-eps(n=10) ~= 0.90
    ## using upper tail:
    1.4826 * qchisq(eps, p, lower.tail=FALSE) / qchisq(0.5, p)
}

## Default wgtFUN() constructors for covComed():
.wgtFUN.covComed <-
    list("01.original" = function(p, ...) {
	     cMah <- .COM.adaptWgt.c(p=p, eps = 0.05)# 1 - eps = 0.95
	     function(d) as.numeric(d < median(d)*cMah) },
	 "01.flex" = function(p, n, control) { ## 'beta' instead of 0.95
	     stopifnot(is.1num(beta <- control$beta), 0 <= beta, beta <= 1)
	     cMah <- 1.4826 * qchisq(beta, p) / qchisq(0.5, p)
	     function(d) as.numeric(d < median(d)*cMah) },
	 "01.adaptive" = function(p, n, ...) { ## 'beta_n' instead of 0.975
	     cMah <- .COM.adaptWgt.c(n,p)
	     function(d) as.numeric(d < cMah) },
	 "sm1.flex" = function(p, n, control) { ## 'beta' / smooth weight
	     stopifnot(is.1num(beta <- control$beta), 0 <= beta, beta <= 1)
	     cMah <- 1.4826 * qchisq(beta, p) / qchisq(0.5, p)
	     function(d) smoothWgt(d / median(d), c=cMah, h = 1) },
	 "sm1.adaptive" = function(p, n, ...) {
	     cMah <- .COM.adaptWgt.c(n=n, p=p)
	     function(d) smoothWgt(d / median(d), c = cMah, h = 1) },
	 "sm2.adaptive" = function(p, n, ...) {
	     cMah <- .COM.adaptWgt.c(n=n, p=p)
	     function(d) smoothWgt(d / median(d), c = cMah, h = 2) }
	 )


comedian <- function (rho, Z, X)
{
    p <- ncol(X)
    U <- svd(rho, nv = 0L)$u
    madX <- apply(X, 2, mad)
    I.mad <- 1/madX
    ## D <- diag(madX)
    ## Q <- D %*% U
    Q <- madX * U
    ## invQ <- solve(Q)
    t.inv.Q <- I.mad * U # = t(solve(Q)) = t(t(U) * D) == t(D) U = D U
    Z <- X %*% t.inv.Q ## for (i in 1:n)  Z[i,] <- invQ %*% X[i,]
    madZ <- apply(Z, 2, mad)
    list(Q=Q, Z=Z, S. = tcrossprod(Q * rep(madZ, each=p))) ## better than
    ##             S. = Q %*% diag(madZ)^2 %*% t(Q)
}

COM <- function(X)
{
    ## Comedian *with* consistency factor.  Falk(1997) was without it.

    stopifnot(is.1num(p <- ncol(X)), p >= 1)
    med <- colMedians(X)
    Y <- sweep(X, 2L, med, `-`)
    COM <- matrix(0., p,p)
    madY <- numeric(p)
    for(i in 1:p) {
        madY[i] <- madYi <- mad(Yi <- Y[,i])
        for(j in seq_len(i-1)) { # j <= i ==> madY[j] "exists"
            COM[j,i] <- COM[i,j] <- median(Yi * Y[,j]) / (madYi * madY[j])
            ## COM[i,j] <- median((Y[,i])*(Y[,j]))
            ## COM[i,j] <- (1.4826^2)*median((Y[,i])*(Y[,j]))
        }
        ## j == i :
        COM[i,i] <- median(Yi^2) / (madYi^2)
    }
    ## return [ 1.4826 = formals(mad)$constant = consistency factor of mad()]
    1.4826^2 * COM
}
