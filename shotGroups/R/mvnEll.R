#####---------------------------------------------------------------------------
## multivariate normal offset ellipse/circle probabilities
#####---------------------------------------------------------------------------

## e characterizes the integration ellipsoid: (x-x0)' e (x-x0) < r^2 -> e = S^{-1}
## cdf
pmvnEll <-
function(r=1, sigma=diag(2), mu, e, x0, lower.tail=TRUE) {
    if(missing(mu)) { mu <- numeric(ncol(sigma)) }
    if(missing(x0)) { x0 <- numeric(ncol(sigma)) }
    if(missing(e))  { e  <- diag(ncol(sigma)) }
    if(!isTRUE(all.equal(as.matrix(sigma), t(sigma)))) { stop("sigma must be symmetric") }
    if(!isTRUE(all.equal(as.matrix(e), t(e))))         { stop("e must be symmetric") }

    ## check e, sigma positive definite
    eEV <- eigen(e)$values
    sEV <- eigen(sigma)$values
    if(!all(eEV >= -sqrt(.Machine$double.eps) * abs(eEV[1]))) {
        stop("e is numerically not positive definite")
    }

    if(!all(sEV >= -sqrt(.Machine$double.eps) * abs(sEV[1]))) {
        stop("sigma is numerically not positive definite")
    }

    ## check dimensions match
    stopifnot(length(x0) == length(mu),
              length(x0) == ncol(e),
              length(x0) == ncol(sigma))

    pp   <- numeric(length(r))               # initialize probabilities to 0
    keep <- which((r > 0) & is.finite(r))    # keep positive r

    ## 1: Mahalanobis transform wrt e -> integrate over unit disc
    L      <- chol(e, pivot=TRUE)
    Esqrt  <- L[ , order(attr(L, "pivot"))]
    xmu1   <- Esqrt %*% (x0-mu)
    sigma1 <- Esqrt %*% sigma %*% t(Esqrt)

    ## 2: rotate with eigenvectors of sigma1 -> decorrelate normal
    S1eig <- eigen(sigma1)
    xmu2  <- t(S1eig$vectors) %*% xmu1

    ## non-centrality parameters
    ncp <- xmu2^2 / S1eig$values
    cqf <- vapply(r[keep], function(x) {
                      CompQuadForm::farebrother(x^2, lambda=S1eig$values, delta=ncp)$res },
                  numeric(1))

    ## NA, NaN, -Inf, Inf (-Inf, Inf will be changed hereafter)
    pp[!is.finite(r)] <- NA_real_

    ## CompQuadForm returns survival probabilities (1-F)
    if(lower.tail) {
        pp[keep] <- 1-cqf
        ## special cases not caught so far
        pp[which(r == -Inf)] <- 0
        pp[which(r ==  Inf)] <- 1
    } else {
        pp[keep] <- cqf
        ## special cases not caught so far
        pp[which(r < 0)]     <- 1
        pp[which(r ==  Inf)] <- 0
    }

    return(pp)
}

## quantile function
qmvnEll <-
function(p, sigma=diag(2), mu, e, x0, lower.tail=TRUE, loUp=NULL) {
    if(missing(mu)) { mu <- numeric(ncol(sigma)) }
    if(missing(x0)) { x0 <- numeric(ncol(sigma)) }
    if(missing(e))  { e  <- diag(ncol(sigma)) }
    if(!isTRUE(all.equal(as.matrix(sigma), t(sigma)))) { stop("sigma must be symmetric") }
    if(!isTRUE(all.equal(as.matrix(e), t(e))))         { stop("e must be symmetric") }

    ## checks on mu, sigma, e are done in getGrubbsParam(), pmvnEll()
    ## initialize quantiles to NA
    qq   <- rep(NA_real_, length(p))
    keep <- which((p >= 0) & (p < 1))
    if(length(keep) < 1) { return(qq) }

    ## determine search interval(s) for uniroot()
    if(is.null(loUp)) {                  # no search interval given
        ## use Grubbs chi^2 quantile for setting root finding interval
        ## Grubbs-Liu chi^2 and actual distribution can diverge
        GP <- getGrubbsParam(sigma=sigma, ctr=(x0-mu), accuracy=TRUE)
        qGrubbs   <- qChisqGrubbs(p[keep], m=GP$m, v=GP$v, muX=GP$muX,
                                  varX=GP$varX, l=GP$l, delta=GP$delta,
                                  lower.tail=lower.tail, type="Liu")
        qGrubbs.6 <- qChisqGrubbs(0.6, m=GP$m, v=GP$v, muX=GP$muX,
                                  varX=GP$varX, l=GP$l, delta=GP$delta,
                                  lower.tail=lower.tail, type="Liu")
        qLo  <- ifelse(p[keep] <= 0.5, 0,         0.25*qGrubbs)
        qUp  <- ifelse(p[keep] <= 0.5, qGrubbs.6, 1.75*qGrubbs)
        loUp <- split(cbind(qLo, qUp), seq_along(p))
    } else {
        if(is.matrix(loUp)) {
            loUp <- split(loUp, seq_len(nrow(loUp)))
        } else if(is.vector(loUp)) {
            loUp <- list(loUp)
        } else if(!is.list(loUp)) {
            stop("loUp must be a list, a matrix, a vector, or missing entirely")
        }
    }

    cdf <- function(r, p, x0, e, mu, sigma, lower.tail) {
        pmvnEll(r=r, sigma=sigma, mu=mu, e=e, x0=x0, lower.tail=lower.tail) - p
    }

    getQ <- function(p, x0, e, mu, sigma, loUp, lower.tail) {
        tryCatch(uniroot(cdf, interval=loUp, p=p, x0=x0, e=e, mu=mu,
                         sigma=sigma, lower.tail=lower.tail)$root,
                 error=function(e) return(NA_real_))
    }

    qq[keep] <- unlist(Map(getQ, p=p[keep], x0=list(x0), e=list(e),
                           mu=list(mu), sigma=list(sigma), loUp=loUp[keep],
                           lower.tail=lower.tail[1]))

    return(qq)
}

## random variates
rmvnEll <-
function(n, sigma=diag(2), mu, e, x0, method=c("eigen", "chol", "cdf"), loUp=NULL) {
    if(missing(mu)) { mu <- numeric(ncol(sigma)) }
    if(missing(x0)) { x0 <- numeric(ncol(sigma)) }
    if(missing(e))  { e  <- diag(ncol(sigma)) }
    if(!isTRUE(all.equal(as.matrix(sigma), t(sigma)))) { stop("sigma must be symmetric") }
    if(!isTRUE(all.equal(as.matrix(e), t(e))))         { stop("e must be symmetric") }

    method <- match.arg(method)

    ## checks on mu, sigma, e are done in getGrubbsParam(), pmvnEll()
    ## if n is a vector, its length determines number of random variates
    n <- if(length(n) > 1L) { length(n) } else { n }

    rn <- if(method == "eigen") {
        lambda <- eigen(sigma)$values    # eigenvalues

        ## simulated 2D normal vectors with mean 0
        X  <- matrix(rnorm(n*ncol(sigma)), nrow=n)     # with identity cov-mat
        xy <- X %*% diag(sqrt(lambda), length(lambda)) # with cov-mat sigma

        ## move to mean mu-x0
        xyMove <- sweep(xy, 2, mu-x0, FUN="+")
        sqrt(rowSums(xyMove^2))          # distances to center
    } else if(method == "chol") {
        CF    <- chol(sigma, pivot=TRUE) # Cholesky-factor
        idx   <- order(attr(CF, "pivot"))
        CFord <- CF[, idx]

        ## simulated 2D normal vectors with mean 0 and cov-mat sigma
        xy <- matrix(rnorm(n*ncol(sigma)), nrow=n) %*% CFord

        ## move to mean mu-x0
        xyMove <- sweep(xy, 2, mu-x0, FUN="+")
        sqrt(rowSums(xyMove^2))          # distances to center
    } else if(method == "cdf") {
        ## root finding of pmvnEll() given uniform random probabilities:
        ## find x such that F(x) - U = 0
        cdf <- function(x, u, sigma, mu, e, x0) {
            pmvnEll(x, sigma=sigma, mu=mu, e=e, x0=x0) - u
        }

        ## find quantile via uniroot() with error handling
        getQ <- function(u, sigma, mu, e, x0, loUp) {
            tryCatch(uniroot(cdf, interval=loUp, u=u, sigma=sigma, mu=mu, e=e, x0=x0)$root,
                     error=function(e) return(NA_real_))
        }

        u <- runif(n)                        # uniform random numbers

        ## determine search interval(s) for uniroot()
        if(is.null(loUp)) {                  # no search interval given
            ## use Grubbs chi^2 quantile for setting root finding interval
            ## Grubbs-Liu chi^2 and actual distribution can diverge
            GP <- getGrubbsParam(sigma=sigma, ctr=(x0-mu), accuracy=TRUE)
            qGrubbs   <- qChisqGrubbs(u, m=GP$m, v=GP$v, muX=GP$muX,
                                      varX=GP$varX, l=GP$l, delta=GP$delta, type="Liu")
            qGrubbs.6 <- qChisqGrubbs(0.6, m=GP$m, v=GP$v, muX=GP$muX,
                                      varX=GP$varX, l=GP$l, delta=GP$delta, type="Liu")
            qLo  <- ifelse(u <= 0.5, 0,         0.25*qGrubbs)
            qUp  <- ifelse(u <= 0.5, qGrubbs.6, 1.75*qGrubbs)
            loUp <- split(cbind(qLo, qUp), seq_along(u))
        } else {
            if(is.matrix(loUp)) {
                loUp <- split(loUp, seq_len(nrow(loUp)))
            } else if(is.vector(loUp)) {
                loUp <- list(loUp)
            } else if(!is.list(loUp)) {
                stop("loUp must be a list, a matrix, a vector, or missing entirely")
            }
        }

        unlist(Map(getQ, u=u, sigma=list(sigma), mu=list(mu), e=list(e),
                   x0=list(x0), loUp=loUp))
    }

    return(rn)
}
