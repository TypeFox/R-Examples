##### need the Matrix package for sparse matrix operation
## require(Matrix)

#### generate data from BMME
#### Input:
####      tgrid: time grid, vector of time points
####      dim  : dimension of Brownian Motion (Wiener process) 
####      sigma: vector of dim, Brownian motion sd, recyclable
####      delta: vector of dim, measurement error sd, recyclable
#### Output:
####      a matrix of dim + 1 columns: time, location coordinates
rbmme <- function(time, dim = 2,  sigma = 1, delta = 1) {  
    n <- length(time)
    dat <- matrix(NA_real_, n, dim)
    t.inc.root <- sqrt(diff(time))
    sigma <- rep(sigma, length = dim)
    delta <- rep(delta, length = dim)
    for (i in 1:dim) {
        bm <- c(0, cumsum(rnorm(n - 1, 0, sd = t.inc.root * sigma[i])))
        err <- rnorm(n, 0, sd = delta[i])
        dat[,i] <- bm + err
    }
    cbind(time = time, dat)
}

#### generate the covariance matrix in sparse structure
#### Input:
####      tinc:  vector of time increment
####      param: vector of (sigma, delta)
#### Output:
####      a banded sparse covariance matrix of the increments
getSparseSigma <- function(tinc, param) {
    n <- length(tinc)
    d0 <- c(tinc * param[1]^2 + 2 * param[2]^2)   
    d1 <- rep(- param[2]^2, n - 1)
    bandSparse(n, k = c(0, 1), diagonals = list(d0, d1), symmetric = TRUE)
}

#### multivariate density using sparse covariance matrix
#### Input:
####      x:     normal vector at which the density needs to be evaluated
####      sigl:  sigl = t(chol(Sigma)), Sigma is sparse covariance matrix
####      log:   if TRUE, return log density
#### Output:
####      log density of density
dmvnormSparse <- function(x, sigl, log = TRUE) {
    n <- NROW(x)
    ## sigl <- t(chol(Sigma))
    part1 <- - n / 2 * log (2 * pi)  - sum(log(diag(sigl)))
    y <- forwardsolve(sigl, x)
    part2 <- - sum(y^2) / 2
    if (log) return(part1 + part2)
    else return(exp(part1 + part2))
}

#### negative loglikelihood which can be fed to optim
#### Input:
####      param: vector of (sigma, delta)
####      dinc:  increment matrix of three columns: time, location_x, location_y
#### Output:
####      loglikelihood
nllk.bmme <- function(param, dinc) {
    if (any(param <= 0)) return(NaN)
    dim <- ncol(dinc) - 1
    Sigma <- getSparseSigma(dinc[,1], c(param[1], param[2]))
    sigl <- t(chol(Sigma))
    llk <- 0
    for (i in 1:dim) {
        llk <- llk + dmvnormSparse(dinc[, i + 1], sigl, log=TRUE)
    }
    return(-llk)

}

#### obtain initial value for sigma and delta by method of moment
#### Input:
####      dat: generated data matrix from gendata ignoring measurement error
####           assuming equal sigma in two directions
#### Output
####      a vector containing rough initial value of sigma 
bmme.start <- function(dat) {
    dif <- apply(dat, 2, diff)
    dim <- ncol(dat) - 1
    st <- sqrt(sum(dif[,-1]^2) / (dim * sum(2 + dif[,1])))
    c(st, st)
}

#### find the MLE by feeding negloglik.full to optim
#### Input:
####      dat:   generated data matrix from gendata
####      start: starting value of (sigma, delta)
#### Output
####      vector containing parameter estimate, standard error, and convergence code

fitBmme <- function(dat, start = NULL, method = "Nelder-Mead",
                    optim.control = list()) {
    if (is.null(start)) start <- bmme.start(dat)
    dinc <- apply(dat, 2, diff)
    fit <- optim(start, nllk.bmme, dinc = dinc, hessian = TRUE, method=method, control = optim.control)
    ## Sigma <- getSparseSigma(dat[,1], fit$par)
    ans <- list(estimate = fit$par,
                var.est = solve(fit$hessian),
                loglik = - fit$value,
                convergence = fit$convergence
                ## vmat.l = t(chol(Sigma))
                )
    ans
}

## remind Vladmir to check linear transformation of (B0, ... Bn, xi0, ..., xin)
getVarObs <- function(tgrid, param) {
    m <- length(tgrid)
    sigma2 <- param[1]^2
    delta2 <- param[2]^2
    diags <- sigma2 * tgrid + delta2
    vmat <- diag(diags)
    offdiags <- sigma2 * rep(tgrid[-m], (m-1):1)
    vmat[row(vmat) > col(vmat)] <- offdiags
    vmat <- vmat + t(vmat) - diag(diags)
    vmat
}
#### density at point x at time t, Prob(X(t) \in dx)
#### Input:
####      times: vector of times at which the density is needed
####      x:     one point at which the density is needed
####      param: vector of (sigma, delta)
####      dat:   data matrix
#### Output:
####      a vector containing the density evaluation at times for point x

dbmme.t.x <- function(times, x, param, SigmaZ.l, dat) {
    dim <- ncol(dat) - 1
    stopifnot(length(x) == dim)
    sigma2 <- param[1]^2
    delta2 <- param[2]^2
    tgrid <- dat[,1]
    sig12 <- sigma2 * outer(times, tgrid, pmin) 
    SigmaZ.l.sig21 <- forwardsolve(SigmaZ.l, t(sig12))
    dens <- 1
    v.t <- sigma2 * times - colSums(SigmaZ.l.sig21^2)
    
    for (i in 1:dim) {
        mu.t <- colSums(SigmaZ.l.sig21 * forwardsolve(SigmaZ.l, dat[, i + 1]))
        dens <- dens * dnorm(x[i], mu.t, sqrt(v.t))
    }
    dens
}

dbb.t.x <- function(times, x, param, dat) {
    dim <- ncol(dat) - 1
    sigma2 <- param[1]^2
    tgrid <- dat[,1]
    t0 <- tgrid[1]
    t1 <- tgrid[2]
    v.t <- (t1 - times) * (times - t0) / (t1 - t0) * sigma2
    dens <- 1
    for (i in 1:dim) {
        mu.t <- dat[1, i + 1] + (times - t0) / (t1 - t0) * diff(dat[, i + 1])
        dens <- dens * dnorm(x[i], mu.t, sqrt(v.t))
    }
    dens
}

## x is a single point
dbmme.x <- function(x, tint, param, SigmaZ.l, dat) {
    tall <- tint[2] - tint[1]
    integrand <- function(u) 100 * dbmme.t.x(u, x, param, SigmaZ.l, dat)
    integrate(integrand, tint[1], tint[2])$value / tall / 100
}

dbmme.x.my <- function(x, tint, param, SigmaZ.l, dat, nstep=100) {
    tall <- tint[2] - tint[1]
    tstep <- tall / nstep
    times <- seq(tint[1]+ tstep / 2, tint[2], by = tstep)
    mean(dbmme.t.x(times, x, param, SigmaZ.l, dat))
}
    
dbb.x <- function(x, tint, param, dat) {
    tall <- tint[2] - tint[1]
    integrand <- function(u) 100 * dbb.t.x(u, x, param, dat)
    integrate(integrand, tint[1], tint[2])$value / tall / 100
}

dbb.x.my <- function(x, tint, param, dat, nstep=100) {
    tall <- tint[2] - tint[1]
    tstep <- tall / nstep
    times <- seq(tint[1]+ tstep / 2, tint[2], by = tstep)
    mean(dbb.t.x(times, x, param, dat))
}

dbmme.x.v <- function(x, tint, param, SigmaZ.l, dat) {
    if (!is.matrix(x)) x <- as.matrix(x)
    apply(x, 1, dbmme.x, tint=tint, param=param, SigmaZ.l=SigmaZ.l, dat=dat)
}

doccuTime <- function(x, tint, param, dat) {
    if (!is.matrix(x)) x <- as.matrix(x)
    SigmaZ <- getVarObs(dat[,1], param)
    SigmaZ.l <- t(chol(SigmaZ))
    doccu <- dbmme.x.v(x, param, SigmaZ.l, dat)
    cbind(x, doccu)
}

