## simulation of breaking time points bbs
## for a 2 state telegraph process
sim1.times.bbz <- function(s, lam1, lam2) {
    tsum <- 0
    tt <- NULL
    while (TRUE) {
        tnew <- rexp(1, lam1) ## starting from state 1
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        if (tsum > s) break
        tnew <- rexp(1, lam2)
        tsum <- tsum + tnew
        tt <- c(tt, tsum)
        if (tsum > s) break
    }
    tt
}

## simulation of a realization given breaking times 
sim1.bbz <- function(s, sigma, time, brtimes, t0moving=TRUE) {
#### time: time points in [0, s]
    tt <- sort(unique(c(time, brtimes)))
    nt <- length(tt)
    nb <- length(brtimes) 
    x <- rep(NA, nt)
    x[1] <- 0
    status <- as.integer(t0moving)
    tend <- brtimes[1]
    j <- 1
    for (i in 2:nt) {
        if (tt[i] <= tend) { ## status unchanged
            if (status == 1)   ## moving
              x[i] <- x[i - 1] + rnorm(1, sd = sigma * sqrt(tt[i] - tt[i - 1]))
            else               ## resting
              x[i] <- x[i - 1]
        }
        if (tt[i] == tend) {
            status <- 1 - status ## switch status
            j <- j + 1
            tend <- brtimes[j]
        }
    }
    x[tt %in% time]
}

## simulation a moving-resting path given a grid time
rMovRes <- function(time, lamM, lamR, sigma, s0, dim = 2) {
    stopifnot(s0 %in% c("m", "r"))
    t0moving <- (s0 == "m")
    lam1 <- if (t0moving) lamM else lamR
    lam2 <- if (t0moving) lamR else lamM
    s <- tail(time, 1)
    brtimes <- sim1.times.bbz(s, lam1, lam2)
    coord <- replicate(dim, sim1.bbz(s, sigma, time, brtimes, t0moving))
    data.frame(time = time, coord)
}


## atom at w = t for dtm.m or dtm.r
dt.atm <- function(t, lam) {
    exp(- lam * t)
}

## dt{mr}.{mr}{mr}
## prob/dens of time spent in state (m or r) by time t,
## the starting state is m or r and ending state is m or r.
##
## dt{mr}.{mr}
## prob/dens of time spent in state (m or r) by time t,
## the starting state is m or r.
##
## using besselI function is much faster than my C code
dtm.mm.b <- function(w, t, lamM, lamR) {
    ## The formula of Zacks (2004, p.500, JAP:497-507) mistakenly has
    ## 2 inside the square root!
    lm <- lamM * w
    lr <- lamR * (t - w)
    elmr <- exp(-lm - lr)
    z <- 2 * sqrt(lm * lr)
    ## ifelse(w > t | w < 0, 0, 
    ##        ifelse(w == t, atm(t, lamM),
    ##               elmr * sqrt(lamM * lamR * w / (t - w)) * besselI(z, 1)))
    ifelse(w > t | w < 0, 0,  elmr * sqrt(lamM * lamR * w / (t - w)) * besselI(z, 1))
}

dtm.mr.b <- function(w, t, lamM, lamR) {
    lm <- lamM * w
    lr <- lamR * (t - w)
    elmr <- exp(-lm - lr)
    z <- 2 * sqrt(lm * lr)
    ifelse(w > t | w < 0, 0, lamM * elmr * besselI(z, 0))
}

## calling C code yields to using modified besselI (mb) as default
dtm.mm <- function(w, t, lamM, lamR, mb = TRUE) { 
    if (mb) return(dtm.mm.b(w, t, lamM, lamR))
    wlen <- length(w)
    t <- rep(t, length.out = wlen)
    dens <- .C("pmm", as.double(w), as.double(t), as.double(lamM), as.double(lamR), as.integer(wlen), dens = double(wlen), PACKAGE="smam")$dens
    ifelse(w > t | w < 0, 0, dens)
    ##           ifelse(w == t, atm(t, lamM), dens))
}

dtm.mr <- function(w, t, lamM, lamR, mb = TRUE) {
    if (mb) return(dtm.mr.b(w, t, lamM, lamR))
    wlen <- length(w)
    t <- rep(t, length.out = wlen)
    dens <- .C("pmr", as.double(w), as.double(t), as.double(lamM), as.double(lamR), as.integer(wlen), dens = double(wlen), PACKAGE="smam")$dens
    ifelse(w > t | w < 0, 0, dens)
}

dtr.rr <- function(w, t, lamM, lamR, mb = TRUE) {
    if (mb) return(dtm.mm.b(w, t, lamR, lamM))
    wlen <- length(w)
    t <- rep(t, length.out=wlen)
    dens <- .C("pmm", as.double(w), as.double(t), as.double(lamR), as.double(lamM), as.integer(wlen), dens = double(wlen), PACKAGE="smam")$dens
    ifelse(w > t | w < 0, 0, dens)
    ##           ifelse(w == t, atm(t, lamR), dens))
}

dtr.rm <- function(w, t, lamM, lamR, mb = TRUE) {
    if (mb) return(dtm.mr.b(w, t, lamR, lamM))
    wlen <- length(w)
    t <- rep(t, length.out=wlen)
    dens <- .C("pmr", as.double(w), as.double(t), as.double(lamR), as.double(lamM), as.integer(wlen), dens = double(wlen), PACKAGE="smam")$dens
    ifelse(w > t | w < 0, 0, dens)
}

dtm.m <- function(w, t, lamM, lamR) {
    dtm.mm(w, t, lamM, lamR) + dtm.mr(w, t, lamM, lamR)
}

dtm.r <- function(w, t, lamM, lamR) {
    dtr.rm(t - w, t, lamM, lamR) + dtr.rr(t - w, t, lamM, lamR)
}
    
dtr.m <- function(w, t, lamM, lamR) {
    dtm.mm(t - w, t, lamM, lamR) + dtm.mr(t - w, t, lamM, lamR)
}

dtr.r <- function(w, t, lamM, lamR) {
    dtr.rm(w, t, lamM, lamR) + dtr.rr(w, t, lamM, lamR)
}

dtm <- function(w, t, lamM, lamR, s0 = NULL) {
    if (is.null(s0)) { ## unconditional
        pm <- lamR / (lamM + lamR)
        pr <- 1 - pm
        pm * dtm.m(w, t, lamM, lamR) + pr * dtm.r(w, t, lamM, lamR)
    }
    else { ## conditional on s0
        if (s0 == "m") dtm.m(w, t, lamM, lamR)
        else dtm.r(w, t, lamM, lamR)
    }
}

dtr <- function(w, t, lamM, lamR, s0 = NULL) {
    if (is.null(s0)) { ## unconditional
        pm <- lamR / (lamM + lamR)
        pr <- 1 - pm
        pm * dtr.m(w, t, lamM, lamR) + pr * dtr.r(w, t, lamM, lamR)
    }
    else { ## conditional on s0
        if (s0 == "m") dtr.m(w, t, lamM, lamR)
        else dtr.r(w, t, lamM, lamR)
    }
}

## density of process at time s, in either state (m or r), given that
## the starting state is m or r

hmm <- function(x, s, lamM, lamR, sigma) {
    ## x is a matrix, s is a match vector
    ## lamM, lamR, and sigma are scalars.
    if (!is.matrix(x)) x <- as.matrix(x) ## make sure x is a matrix
    d <- ncol(x) ## can be univariate or bivariate
    n <- nrow(x)
    s <- rep(s, length.out=n)
    term1 <- exp(-lamM * s) * apply(dnorm(x, sd = sigma * sqrt(s)), 1, prod)
    term2 <- sapply(1:n,
                    function(i) {
                        f <- function(w) {
                            t1 <- prod(dnorm(x[i,,drop=FALSE], sd=sigma*sqrt(w)))
                            t2 <- dtm.mm(w, s[i], lamM, lamR)
                            t1 * t2
                        }
                        vf <- Vectorize(f)
                        integrate(vf, 0, s[i])$value
                    })
    term1 + term2
}

hmr <- function(x, s, lamM, lamR, sigma) {
    ## x is a matrix, s is a match vector
    ## lamM, lamR, and sigma are scalars.
    if (!is.matrix(x)) x <- as.matrix(x) ## make sure x is a matrix
    d <- ncol(x) ## can be univariate or bivariate
    n <- nrow(x)
    s <- rep(s, length.out=n)
    sapply(1:n,
           function(i) {
               f <- function(w) {
                   t1 <- prod(dnorm(x[i,,drop=FALSE], sd=sigma*sqrt(w)))
                   t2 <- dtm.mr(w, s[i], lamM, lamR)
                   t1 * t2
               }
               vf <- Vectorize(f)
               integrate(vf, 0, s[i])$value
           })
}


hrr <- function(x, s, lamM, lamR, sigma) {
    ## x is a matrix, s is a match vector
    ## lamM, lamR, and sigma are scalars.
    if (!is.matrix(x)) x <- as.matrix(x) ## make sure x is a matrix
    d <- ncol(x) ## can be univariate or bivariate
    n <- nrow(x)
    s <- rep(s, length.out=n)
    sapply(1:n,
           function(i) {
               f <- function(w) {
                   t1 <- prod(dnorm(x[i,,drop=FALSE], sd=sigma*sqrt(s[i] - w)))
                   t2 <- dtr.rr(w, s[i], lamM, lamR)
                   t1 * t2
               }
               vf <- Vectorize(f)
               integrate(vf, 0, s[i])$value
           })
}

hrm <- function(x, s, lamM, lamR, sigma) {
    ## x is a matrix, s is a match vector
    ## lamM, lamR, and sigma are scalars.
    if (!is.matrix(x)) x <- as.matrix(x) ## make sure x is a matrix
    d <- ncol(x) ## can be univariate or bivariate
    n <- nrow(x)
    s <- rep(s, length.out=n)
    sapply(1:n,
           function(i) {
               f <- function(w) {
                   t1 <- prod(dnorm(x[i,,drop=FALSE], sd=sigma*sqrt(s[i] - w)))
                   t2 <- dtr.rm(w, s[i], lamM, lamR)
                   t1 * t2
               }
               vf <- Vectorize(f)
               integrate(vf, 0, s[i])$value
           })
}

hm <- function(x, s, lamM, lamR, sigma) {
    ## integrate(hm, -infty, +infty) = 1 
    hmm(x, s, lamM, lamR, sigma) + hmr(x, s, lamM, lamR, sigma)
}

hr <- function(x, s, lamM, lamR, sigma) {
    ## integrate(hr, -infty, +infty) = 1 - exp(-lamR * s)
    hrm(x, s, lamM, lamR, sigma) + hrr(x, s, lamM, lamR, sigma)
}

## marginal likelihood for inferences
## marginal density of one increment
## x: the increment (not the original locations)
dmrbm.m1 <- function(x, tt, lamM, lamR, sigma) {
    if (lamM <= 0 || lamR <=0 || sigma <= 0) return(NA)
    if (!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    ## delta <- apply(x, 1, prod) == 0
    delta <- apply(x == 0, 1, all) ## staying
    pm <- 1 / lamM / (1 / lamM + 1 / lamR)
    pr <- 1 - pm
    tt <- rep(tt, length.out = n)
    dens <- pr * exp(-lamR * tt)
    if (!all(delta)) {
        x <- x[!delta, , drop=FALSE]
        tt <- tt[!delta]
        hm <- hmm(x, tt, lamM, lamR, sigma) + hmr(x, tt, lamM, lamR, sigma)
        hr <- hrm(x, tt, lamM, lamR, sigma) + hrr(x, tt, lamM, lamR, sigma)
        dens[!delta] <- pm * hm + pr * hr
    }
    dens
}

## likelihood for given one set of parameter value
cllk.m1.p <- function(x, t, lamM, lamR, sigma) { 
    cllk <- try(sum(log(dmrbm.m1(x, t, lamM, lamR, sigma))), silent=TRUE)
    if (inherits(cllk, "try-error")) NA else cllk
}


cllk.m1.p2v <- Vectorize(cllk.m1.p, c("lamM", "lamR", "sigma"))

cllk.m1 <- function(theta, data) {
    lamM <- theta[1]
    lamR <- theta[2]
    sigma <- theta[3]
    dinc <- apply(data, 2, diff)
    cllk.m1.p2v(dinc[,-1], dinc[,1], lamM, lamR, sigma)    
}

ncllk.m1.inc <- function(theta, data) { ## data is increment already
    lamM <- theta[1]
    lamR <- theta[2]
    sigma <- theta[3]
    -cllk.m1.p2v(data[,-1], data[,1], lamM, lamR, sigma)    
}


fitMovRes <- function(data, start, method = "Nelder-Mead",
                      optim.control = list()) {
    dinc <- apply(data, 2, diff)
    fit <- optim(start, ncllk.m1.inc, data = dinc,
                 method = method,
                 control = optim.control)
    list(estimate = fit$par,
         cloglik = -fit$value,
         convergence = fit$convergence)
}

## dx <- function(x, t, lamM, lamR, Sigma) {
##     pm <- 1 / lamM / (1 / lamM + 1 / lamR)
##     pr <- 1 - pm
##     hm <- dxmm(x, t, lamM, lamR, Sigma) + dxmr(x, t, lamM, lamR, Sigma)
##     hr <- dxrm(x, t, lamM, lamR, Sigma) + dxrr(x, t, lamM, lamR, Sigma)
##     pm * hm + pr * hr
## }

## cllk.1inc <- function(theta, xinc, t) {
##     Sigma <- diag(rep(theta[3], 2))
##     sum(log(dx(xinc, t, theta[1], theta[2], Sigma)))
## }


## #################################################################
## ## old code
## #################################################################
## dmnorm <- function(x, Sigma, w) {
##   ## x is a single observation vector
##   ## w is a vector
##   ## need to return a vector of the same length as w
##   if (!is.matrix(x)) x <- as.matrix(x)
##   if (!is.matrix(Sigma)) Sigma <- as.matrix(Sigma)
##   sapply(w, function(t) dmvnorm(x, sigma = t * Sigma))
## }
   
## dxmm <- function(x, s, lamM, lamR, Sigma) {
##   if (is.vector(x)) x <- matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   s <- rep(s, length.out=xlen)
##   val <- rep(NA, xlen)
##   for (i in 1:xlen) {
##     f <- function(w) dmnorm(matrix(x[i, ], 1, d), Sigma, w) * dtm.mm(w, s[i], lamM, lamR)
##     val[i] <- integrate(f, 0, s[i])$value
##     val[i] <- val[i] + exp(-lamM * s[i]) * dmnorm(matrix(x[i, ], 1, d), Sigma, s[i])
##   }
##   val
## }

## ## hmr(x, t)
## dxmr <- function(x, s, lamM, lamR, Sigma) {
##   if(is.vector(x)) x <- matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   s <- rep(s, length.out=xlen)
##   val <- rep(NA, xlen)
##   for (i in 1:xlen) {
##     f <- function(w) dmnorm(matrix(x[i, ], 1, d), Sigma, w) * dtm.mr(w, s[i], lamM, lamR)
##     val[i] <- integrate(f, 0, s[i])$value
##   }
    
##   val
## }

## ## hrr(x, t)
## dxrr <- function(x, s, lamM, lamR, Sigma) {
##   if(is.vector(x)) x <- matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   s <- rep(s, length.out=xlen)
##   val <- rep(NA, xlen)
##   for (i in 1:xlen) {
##     f <- function(w) dmnorm(matrix(x[i, ], 1, d), Sigma, s[i] - w) * dtr.rr(w, s[i], lamM, lamR)
##     val[i] <- integrate(f, 0, s[i])$value
##   }

##   val
## }

## ## hrm(x, t)
## dxrm <- function(x, s, lamM, lamR, Sigma) {
##   if(is.vector(x)) x <- matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   s <- rep(s, length.out=xlen)
##   val <- rep(NA, xlen)
##   for (i in 1:xlen) {
##     f <- function(w) dmnorm(matrix(x[i, ], 1, d), Sigma, s[i] - w) * dtr.rm(w, s[i], lamM, lamR)
##     val[i] <- integrate(f, 0, s[i])$value
##   }

##   val
## }


## ## conditional density
## ## Pr[ X_s in (x, x + dx), X_t - X_s in (y, y + dy) | s(0) = 1]
## dxxm <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- 1 / lamM / (1 / lamM + 1 / lamR)
##   pr <- 1 - pm
##   dxxm.m <- dxmm(x, s, lamM, lamR, Sigma) * (dxmm(y, t - s, lamM, lamR, Sigma) +
##                                             dxmr(y, t - s, lamM, lamR, Sigma))
##   dxxm.r <- dxmr(x, s, lamM, lamR, Sigma) * (dxrr(y, t - s, lamM, lamR, Sigma) +
##                                             dxrm(y, t - s, lamM, lamR, Sigma))
##   dxxm.m+dxxm.r 
## }

## ## Pr[ X_s in (x, x + dx), X_t - X_s in (y, y + dy) | s(0) = 0]
## dxxr <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- 1 / lamM / (1 / lamM + 1 / lamR)
##   pr <- 1 - pm
##   dxxr.m <- dxrm(x, s, lamM, lamR, Sigma) * (dxmm(y, t - s, lamM, lamR, Sigma) +
##                                             dxmr(y, t - s, lamM, lamR, Sigma))
##   dxxr.r <- dxrr(x, s, lamM, lamR, Sigma) * (dxrr(y, t - s, lamM, lamR, Sigma) +
##                                             dxrm(y, t - s, lamM, lamR, Sigma))
##   dxxr.m+dxxr.r  
## }

## ## Pr[ X_s in (x, x + dx), X_t - X_s in (y, y + dy) | s(0) = 1 or 0]
## dxx <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- 1 / lamM / (1 / lamM + 1 / lamR)
##   pr <- 1 - pm
##   dxxm(x, y, s, t, lamM, lamR, Sigma) * pm + dxxr(x, y, s, t, lamM, lamR, Sigma) * pr  
## }



## BDd <- function(y, t, lamM, lamR, Sigma) {
##   pm <- lamR / (lamR + lamM)
##   pr <- lamM / (lamR + lamM)
##   hm <- dxmm(y, t, lamM, lamR, Sigma) + dxmr(y, t, lamM, lamR, Sigma)
##   hr <- dxrm(y, t, lamM, lamR, Sigma) + dxrr(y, t, lamM, lamR, Sigma)
##   bdd <- pm * hm + pr * hr
##   bdd
## }

## ## pr[ X(u) in dx | X(t) = y] 
## BD <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- lamR / (lamR + lamM)
##   pr <- lamM / (lamR + lamM)

##   if(is.vector(x)) x <- as.matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   zero <- ifelse(d == 1, 0, matrix(0, 1, d))

##   bd <- matrix(NA, length(s), xlen)
##   bdd <- BDd(y, t, lamM, lamR, Sigma)

##   for(j in 1:xlen) {  
    
##      # x != 0, y != 0, x != Y
##     if((!identical(matrix(as.numeric(x[j,]), 1, d), zero)) & (!identical(y, zero)) & (!identical(matrix(as.numeric(x[j,]), 1, d), y))) {
##       bdn <- sapply(1:length(s), function(i) dxx(matrix(x[j,], 1, d), y - matrix(x[j,], 1, d), s[i], t, lamM, lamR, Sigma))
##      # bdn <- sapply(1:length(s), function(i) pm * dxxm(x = matrix(x[j,], 1, d), y = y - matrix(x[j,], 1, d), s[i], t, lamM, lamR, Sigma) + pr * dxxr(x = matrix(x[j,], 1, d), y = y - x[j,], s[i], t, lamM, lamR, Sigma))
##       bd[, j] <- sapply(1:length(s), function(i) ifelse(bdd > 0, bdn[i] / bdd, 0))
##     }
##     # x = 0, y != 0
##     else if ((identical(matrix(as.numeric(x[j,]), 1, d), zero)) & (!identical(y, zero))) {  
##       bd[, j] <- sapply(1:length(s), function(i) ifelse(bdd > 0, (pr * exp(-lamR*(s[i])) * (dxrm(y, t-s[i], lamM, lamR, Sigma)+dxrr(y, t-s[i], lamM, lamR, Sigma))) / bdd, 0))
##     }
##     # x = y, y != 0
##     else if ((identical(matrix(as.numeric(x[j,]), 1, d), y)) & (!identical(y, zero))) {  
##       bd[, j] <- sapply(1:length(s), function(i) ifelse(bdd > 0, (pr * exp(-lamR*(t-s[i])) * (dxrm(y, s[i], lamM, lamR, Sigma)+dxrr(y, s[i], lamM, lamR, Sigma))) / bdd, 0))
##     }
##     # y = 0
##     else if (identical(y, zero)) {
##       bd[, j] <- sapply(1:length(s), function(i) ifelse(identical(matrix(as.numeric(x[j,]), 1, d), zero), 1, 0))
##     }


##   }

##   bd
## }


## BDIni <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- lamR / (lamR + lamM)
##   pr <- lamM / (lamR + lamM)

##   if(is.vector(x)) x <- as.matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   zero <- matrix(0, 1, d)

##   bd <- matrix(NA, length(s), xlen)
##   bdd <- BDd(y, t, lamM, lamR, Sigma)

##   for(j in 1:xlen) {
##     # x = 0, y != 0
##     bd[, j] <- sapply(1:length(s), function(i) ifelse(bdd > 0, (pr * exp(-lamR*(s[i])) * (dxrm(y, t-s[i], lamM, lamR, Sigma)+dxrr(y, t-s[i], lamM, lamR, Sigma))) / bdd, 0))

##   }

##   bd
## }

## BDEnd <- function(x, y, s, t, lamM, lamR, Sigma) {
##   pm <- lamR / (lamR + lamM)
##   pr <- lamM / (lamR + lamM)

##   if(is.vector(x)) x <- as.matrix(x)
##   d <- ncol(x)
##   xlen <- nrow(x)
##   zero <- matrix(0, 1, d)

##   bd <- matrix(NA, length(s), xlen)
##   bdd <- BDd(y, t, lamM, lamR, Sigma)

##   for(j in 1:xlen) {
##     # x = y, y != 0
##     bd[, j] <- sapply(1:length(s), function(i) ifelse(bdd > 0, (pr * exp(-lamR*(t-s[i])) * (dxrm(y, s[i], lamM, lamR, Sigma)+dxrr(y, s[i], lamM, lamR, Sigma))) / bdd, 0))
  
##   }

##   bd
## }
















