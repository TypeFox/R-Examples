
###### Function to create block diagonal matrices (from R-help, I suppose)
bdiag <- function(...)
{
    if (nargs() == 1)
        x <- as.list(...)
    else
        x <- list(...)
    n <- length(x)
    if(n==0) return(NULL)
    x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
                stop("Zero-length component in x"))
    d <- array(unlist(lapply(x, dim)), c(2, n))
    rr <- d[1,]
    cc <- d[2,]
    rsum <- sum(rr)
    csum <- sum(cc)
    out <- array(0, c(rsum, csum))
    ind <- array(0, c(4, n))
    rcum <- cumsum(rr)
    ccum <- cumsum(cc)
    ind[1,-1] <- rcum[-n]
    ind[2,] <- rcum
    ind[3,-1] <- ccum[-n]
    ind[4,] <- ccum
    imat <- array(1:(rsum * csum), c(rsum, csum))
    iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
                                                 (y[3]+1):y[4]], imat=imat)
    iuse <- as.vector(unlist(iuse))
    out[iuse] <- unlist(x)
    return(out)
}

###### From R^p to the AR parameters of a stationary AR(p) (univariate)
###### For the multivariate analog, see Ansley and Kohn, 1986,
###### J. Statist. Comput. Simul. 24.
ARtransPars <- function(raw) {
    p <- length(raw)
    return(.Call("ARtranspar", p, as.double(raw), PACKAGE="dlm"))
}

###### Constructor for dlm objects
dlm <- function(...) {
    if (nargs() == 1)
        x <- as.list(...)
    else
        x <- list(...)
    ## required components
    nm <- c("m0", "C0", "FF", "V", "GG", "W")
    nmInd <- match(nm, names(x))
    if (any(is.na(nmInd)))
        stop(paste("Component(s)", paste(nm[is.na(nmInd)], collapse=", "),
                   "is (are) missing"))
    x[nmInd[-1]] <- lapply(x[nmInd[-1]], as.matrix)
    if (!is.numeric(x$FF))
        stop("Component FF must be numeric")
    m <- nrow(x$FF)
    p <- ncol(x$FF)
    if (!is.numeric(x$V))
        stop("Component V must be numeric")
    if (!( nrow(x$V) == m && ncol(x$V) == m))
        stop("Incompatible dimensions of matrices")
    if (!is.numeric(x$GG))
        stop("Component GG must be numeric")
    if (!( nrow(x$GG) == p && ncol(x$GG) == p))
        stop("Incompatible dimensions of matrices")
    if (!is.numeric(x$W))
        stop("Component W must be numeric")
    if (!( nrow(x$W) == p && ncol(x$W) == p))
        stop("Incompatible dimensions of matrices")
    if (!is.numeric(x$C0))
        stop("Component C0 must be numeric")
    if (!( nrow(x$C0) == p && ncol(x$C0) == p))
        stop("Incompatible dimensions of matrices")
    if (!( is.numeric(x$m0) && NCOL(x$m0) == 1 && NROW(x$m0) == p ))
        stop(paste("Component m0 must be a numeric vector of length",
                   "\n equal to ncol of component FF, or a matrix with one column and",
                   "\n number of rows equal to ncol of component FF"))
    if (!(all.equal(x$C0, t(x$C0)) && all(eigen(x$C0)$values >= 0)))
        stop("C0 is not a valid variance matrix")
    if (any( c(is.na(x$m0), is.na(x$C0))))
        stop("Missing values are not allowed in components m0 and C0")
    ## extra components for time-varying dlm
    nm1 <- c("JFF", "JV", "JGG", "JW")
    nm1Ind <- match(nm1, names(x))
    if (all(is.na(nm1Ind))) {
        if (!(all.equal(x$V, t(x$V)) && all(eigen(x$V)$values >= 0)))
            stop("V is not a valid variance matrix")
        if (!(all.equal(x$W, t(x$W)) && all(eigen(x$W)$values >= 0)))
            stop("W is not a valid variance matrix")
        mod <- x[nmInd]
        class(mod) <- "dlm"
        return(mod)
    }
    x[nm1Ind[!is.na(nm1Ind)]] <-
        lapply(x[nm1Ind[!is.na(nm1Ind)]], function(x) if (!is.null(x)) as.matrix(x))
    if (!is.null(x$JFF)) {
        if (!(is.numeric(x$JFF) &&
              nrow(x$JFF) == m && ncol(x$JFF) == p))
            stop("Invalid component JFF")
        JFF <- round(x$JFF)
        if (all(JFF == 0)) JFF <- NULL
    } else
    JFF <- NULL
    if (!is.null(x$JV)) {
        if (!(is.numeric(x$JV) &&
              nrow(x$JV) == m && ncol(x$JV) == m))
            stop("Invalid component JV")
        JV <- round(x$JV)
        if (all(JV == 0)) JV <- NULL
    } else
    JV <- NULL
    if (!is.null(x$JGG)) {
        if (!(is.numeric(x$JGG) &&
              nrow(x$JGG) == p && ncol(x$JGG) == p))
            stop("Invalid component JGG")
        JGG <- round(x$JGG)
        if (all(JGG == 0)) JGG <- NULL
    } else
    JGG <- NULL
    if (!is.null(x$JW)) {
        if (!(is.numeric(x$JW) &&
              nrow(x$JW) == p && ncol(x$JW) == p))
            stop("Invalid component JW")
        JW <- round(x$JW)
        if (all(JW == 0)) JW <- NULL
    } else
    JW <- NULL
    mx <- max(c(JFF, JV, JGG, JW))
    if ( mx <= 0 ) {
        mod <- x[nmInd]
        class(mod) <- "dlm"
        return(mod)
    }
    if ( is.null(x$X) )
        stop("Component X must be provided for time-varying models")
    x$X <- as.matrix(x$X)
    if (!(is.numeric(x$X) && ncol(x$X) >= mx))
        stop("Invalid component X")
    mod <- c(x[nmInd], list(JFF=JFF, JV=JV, JGG=JGG, JW=JW, X=x$X))
    class(mod) <- "dlm"
    return(mod)
}

is.dlm <- function(obj) inherits(obj, "dlm")

as.dlm <- function(obj)
    if (is.dlm(obj)) obj else dlm(obj)

is.dlmFiltered <- function(obj) inherits(obj, "dlmFiltered")


###### Extractors and replacement functions
FF <- function(x) UseMethod("FF")
"FF<-" <- function(x, value) UseMethod("FF<-")
GG <- function(x) UseMethod("GG")
"GG<-" <- function(x, value) UseMethod("GG<-")
V <- function(x) UseMethod("V")
"V<-" <- function(x, value) UseMethod("V<-")
W <- function(x) UseMethod("W")
"W<-" <- function(x, value) UseMethod("W<-")
C0 <- function(x) UseMethod("C0")
"C0<-" <- function(x, value) UseMethod("C0<-")
m0 <- function(x) UseMethod("m0")
"m0<-" <- function(x, value) UseMethod("m0<-")

FF.dlm <- function(x)
{
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    if (!is.null(x$JFF))
        warning(paste("Time varying", sQuote("F")))
    return(x$FF)
}

"FF<-.dlm" <- function(x, value)
{
    value <- as.matrix(value)
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    if (any(!is.finite(value)))
        stop(paste("finite values needed in", sQuote("F")))
    if (!(is.numeric(value) && is.matrix(value)))
        stop(paste(deparse(substitute(value)), "is not a valid observation matrix"))
    if ( any(dim(x$FF) != dim(value)) )
        stop(paste("wrong dimension of", deparse(substitute(value))))
    if (!is.null(x$JFF))
        warning(paste("time varying", sQuote("F"), "in", sQuote("x")))
    x$FF <- value
    return(x)
}

GG.dlm <- function(x)
{
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    if (!is.null(x$JGG))
        warning(paste("Time varying", sQuote("G")))
    return(x$GG)
}

"GG<-.dlm" <- function(x, value)
{
    value <- as.matrix(value)
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    if (any(!is.finite(value)))
        stop(paste("finite values needed in", sQuote("G")))
    if (!(is.numeric(value) && is.matrix(value)))
        stop(paste(deparse(substitute(value)), "is not a valid evolution matrix"))
    if ( any(dim(x$GG) != dim(value)) )
        stop(paste("wrong dimension of", deparse(substitute(value))))
    if (!is.null(x$JGG))
        warning(paste("Time varying", sQuote("G"), "in", sQuote("x")))
    x$GG <- value
    return(x)
}

V.dlm <- function(x)
{
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    if (!is.null(x$JV))
        warning(paste("Time varying", sQuote("V")))
    return(x$V)
}

"V<-.dlm" <- function(x, value)
{
    value <- as.matrix(value)
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    if (any(!is.finite(value)))
        stop(paste("finite values needed in", sQuote("V")))
    if (!(is.numeric(value) && is.matrix(value) && all.equal(value, t(value))
          && all(eigen(value)$values > -.Machine$double.neg.eps)))
        stop(paste(deparse(substitute(value)), "is not a valid variance matrix"))
    if ( any(dim(x$V) != dim(value)) )
        stop(paste("wrong dimension of", deparse(substitute(value))))
    if (!is.null(x$JV))
        warning(paste("Time varying", sQuote("V"), "in", sQuote("x")))
    x$V <- value
    return(x)
}

W.dlm <- function(x)
{
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    if (!is.null(x$JW))
        warning(paste("Time varying", sQuote("W")))
    return(x$W)
}

"W<-.dlm" <- function(x, value)
{
    value <- as.matrix(value)
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    if (any(!is.finite(value)))
        stop(paste("finite values needed in", sQuote("W")))
    if (!(is.numeric(value) && is.matrix(value) && all.equal(value, t(value))
          && all(eigen(value)$values > -.Machine$double.neg.eps)))
        stop(paste(deparse(substitute(value)), "is not a valid variance matrix"))
    if ( any(dim(x$W) != dim(value)) )
        stop(paste("wrong dimension of", deparse(substitute(value))))
    if (!is.null(x$JW))
        warning(paste("Time varying", sQuote("W"), "in", sQuote("x")))
    x$W <- value
    return(x)
}

C0.dlm <- function(x)
{
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    return(x$C0)
}

"C0<-.dlm" <- function(x, value)
{
    value <- as.matrix(value)
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    if (any(!is.finite(value)))
        stop(paste("finite values needed in", sQuote("C0")))
    if (!(is.numeric(value) && is.matrix(value) && all.equal(value, t(value))
          && all(eigen(value)$values > 0)))
        stop(paste(deparse(substitute(value)), "is not a valid variance matrix"))
    if ( any(dim(x$C0) != dim(value)) )
        stop(paste("wrong dimension of", deparse(substitute(value))))
    x$C0 <- value
    return(x)
}

m0.dlm <- function(x) {
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    return(x$m0)
}

"m0<-.dlm" <- function(x, value)
{
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    if (any(!is.finite(value)))
        stop(paste("finite values needed in", sQuote("m0")))
    if (!is.numeric(value))
        stop(paste(deparse(substitute(value)), "is not numeric"))
    if (!(length(value) == length(x$m0) && NCOL(value) == 1))
        stop(paste("wrong length/dimension of", deparse(substitute(value))))
    x$m0 <- value
    return(x)
}

### Stuff for time-varying models
JFF <- function(x) UseMethod("JFF")
"JFF<-" <- function(x, value) UseMethod("JFF<-")
JGG <- function(x) UseMethod("JGG")
"JGG<-" <- function(x, value) UseMethod("JGG<-")
JV <- function(x) UseMethod("JV")
"JV<-" <- function(x, value) UseMethod("JV<-")
JW <- function(x) UseMethod("JW")
"JW<-" <- function(x, value) UseMethod("JW<-")
X <- function(x) UseMethod("X")
"X<-" <- function(x, value) UseMethod("X<-")

JFF.dlm <- function(x) {
    return(x$JFF)
}

"JFF<-.dlm" <- function(x, value)
{
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    value <- as.matrix(value)
    dm <- dim(value)
    value <- as.integer(value)
    dim(value) <- dm
    if (any(is.na(value)))
        stop(paste("missing values not allowed in", sQuote("JFF")))
    if ( any(dim(x$FF) != dim(value)) )
        stop(paste("wrong dimension of", deparse(substitute(value))))
    x$JFF <- value
    return(x)
}

JGG.dlm <- function(x) {
    return(x$JGG)
}

"JGG<-.dlm" <- function(x, value)
{
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    value <- as.matrix(value)
    dm <- dim(value)
    value <- as.integer(value)
    dim(value) <- dm
    if (any(is.na(value)))
        stop(paste("missing values not allowed in", sQuote("JGG")))
    if ( any(dim(x$GG) != dim(value)) )
        stop(paste("wrong dimension of", deparse(substitute(value))))
    x$JGG <- value
    return(x)
}

JV.dlm <- function(x) {
    return(x$JV)
}

"JV<-.dlm" <- function(x, value)
{
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    value <- as.matrix(value)
    dm <- dim(value)
    value <- as.integer(value)
    dim(value) <- dm
    if (any(is.na(value)))
        stop(paste("missing values not allowed in", sQuote("JV")))
    if ( any(dim(x$V) != dim(value)) )
        stop(paste("wrong dimension of", deparse(substitute(value))))
    x$JV <- value
    return(x)
}

JW.dlm <- function(x) {
    return(x$JW)
}

"JW<-.dlm" <- function(x, value)
{
    if (!is.dlm(x))
        stop(paste(deparse(substitute(x)), "is not a dlm object"))
    value <- as.matrix(value)
    dm <- dim(value)
    value <- as.integer(value)
    dim(value) <- dm
    if (any(is.na(value)))
        stop(paste("missing values not allowed in", sQuote("JW")))
    if ( any(dim(x$W) != dim(value)) )
        stop(paste("wrong dimension of", deparse(substitute(value))))
    x$JW <- value
    return(x)
}

X.dlm <- function(x) {
    return(x$X)
}

"X<-.dlm" <- function(x, value)
{
    value <- as.matrix(value)
    dm <- dim(value)
    value <- as.numeric(value)
    dim(value) <- dm
    if (any(!is.finite(value)))
        stop(paste("finite values needed in", sQuote("X")))
    x$X <- value
    return(x)
}


###### Regression
dlmModReg <- function(X, addInt=TRUE, dV=1, dW=rep(0,NCOL(X)+addInt),
                      m0=rep(0,length(dW)), C0=1e7*diag(nrow=length(dW)))
{
    ## sanity checks
    p <- NCOL(X) + addInt
    if (!( length(dV)==1 && length(dW)==p &&
          length(m0)==p && nrow(C0)==p && ncol(C0)==p))
        stop("Inconsistent dimensions of arguments")
    X <- as.matrix(X)
    JFF <- matrix(1:ncol(X),nrow=1)
    if (addInt)
        JFF <- cbind(0,JFF)
    mod <- list(
                m0 = m0,
                C0 = C0,
                FF = matrix(1,1,p),
                V = as.matrix(dV),
                GG = diag(nrow=p),
                W = diag(x=dW, nrow=p),
                JFF = JFF,
                JV = NULL,
                JGG = NULL,
                JW = NULL,
                X = X)
    class(mod) <- "dlm"
    return(mod)
}

###### Polynomial trend
dlmModPoly <- function(order=2, dV=1,
                       dW=c(rep(0,order-1),1),
                       m0=rep(0,order), C0=1e7*diag(nrow=order))
{
    C0 <- as.matrix(C0)
    ## sanity checks
    if (!( length(dV)==1 && length(dW)==order &&
          length(m0)==order && nrow(C0)==order && ncol(C0)==order))
        stop("Inconsistent dimensions of arguments")
    GG <- diag(order)
    GG[row(GG) == col(GG)-1] <- 1
    mod <- list(
                m0 = m0,
                C0 = C0,
                FF = matrix(c(1,rep(0,order-1)),nrow=1),
                V = as.matrix(dV),
                GG = GG,
                W = diag(dW,nrow=length(dW)),
                JFF = NULL,
                JV = NULL,
                JGG = NULL,
                JW = NULL)
    class(mod) <- "dlm"
    return(mod)
}

###### Seasonal factors
dlmModSeas <- function(frequency, dV=1, dW=c(1,rep(0,frequency-2)),
                       m0=rep(0,frequency-1),
                       C0=1e7*diag(nrow=frequency-1))
{
    frequency <- as.integer(frequency)
    p <- frequency - 1
    if (!( length(dV) == 1 && length(dW) == p && length(m0) == p &&
          nrow(C0) == p && ncol(C0) == p ))
        stop("Inconsistent dimensions of arguments")
    GG <- matrix(0,p,p)
    GG[row(GG) == col(GG)+1] <- 1
    GG[1,] <- -1
    mod <- list(
                m0 = m0,
                C0 = C0,
                FF = matrix(c(1,rep(0,p-1)),nrow=1),
                V = as.matrix(dV),
                GG = GG,
                W = diag(dW,nrow=length(dW)),
                JFF = NULL,
                JV = NULL,
                JGG = NULL,
                JW = NULL)
    class(mod) <- "dlm"
    return(mod)
}

###### Fourier representation
dlmModTrig <- function(s, q, om, tau, dV=1, dW=0, m0, C0)
{
    if ( hasArg(s) ) {
        if ( hasArg(om) )
            stop("Cannot specify both 's' and 'om'")
        if ( hasArg(tau) )
            stop("Cannot specify both 's' and 'tau'")
        s <- as.integer(s)
        sHalf <- s %/% 2
        if ( hasArg(q) ) {
            q <- as.integer(q)
            if ( q > sHalf )
                stop(paste("Can use",sHalf,"components at most (",q,"were asked for )"))
            if ( q <= 0 ) stop("'q' must be positive")
        } else {
            q <- sHalf }
        om <- 2 * base::pi / s
        evenAll <- (q == sHalf) && !(s %% 2)
        if ( sHalf == 1 && evenAll ) {
            FF <- matrix(1,1,1)
            GG <- matrix(-1,1,1)
        } else {
            H1 <- diag(x=cos(om), nrow=2)
            H1[1,2] <- sin(om)
            H1[2,1] <- - H1[1,2]
            if ( q > 1 ) {
                h <- vector("list",q)
                h[[1]] <- H1
                i <- 2
                while ( i < q ) {
                    h[[i]] <- h[[i-1]] %*% H1
                    i <- i + 1
                }
                h[[q]] <- if ( evenAll ) matrix(-1,1,1) else h[[q-1]] %*% H1
                GG <- bdiag(h)
            } else {
                GG <- H1 }
        }
    } else {
        if ( !hasArg(om) )
            if ( !hasArg(tau) )
                stop("One of 's', 'om' or 'tau' must be specified")
            else om <- 2 * base::pi / tau
        if ( !hasArg(q) )
            stop("When 'om' or 'tau' is specified, 'q' must be supplied as well")
        q <- as.integer(q)
        if ( q <= 0 ) stop("'q' must be positive")
        evenAll <- FALSE
        H1 <- diag(x=cos(om), nrow=2)
        H1[1,2] <- sin(om)
        H1[2,1] <- - H1[1,2]
        if ( q > 1 ) {
            h <- vector("list",q)
            h[[1]] <- H1
            for ( i in 2:q )
                h[[i]] <- h[[i-1]] %*% H1
            GG <- bdiag(h)
        } else {
            GG <- H1 }
    }
    if ( hasArg(m0) ) {
        if ( length(m0) != nrow(GG) ) stop("Wrong length of 'm0'")
    } else {
        m0 <- rep(0,nrow(GG)) }
    if ( hasArg(C0) ) {
        if ( (nrow(C0) != nrow(GG)) || (ncol(C0) != nrow(GG)) )
            stop("Wrong dimension of 'C0'")
    } else {
        C0 <- diag(x=1e7, nrow=nrow(GG)) }
    mod <- list(
                m0 = m0,
                C0 = C0,
                FF = matrix(if (evenAll) c(rep(c(1,0),q-1),1)
                else c(1,0),1,nrow(GG)),
                V = matrix(dV,1,1),
                GG = GG,
                W = diag(x=dW,nrow=nrow(GG)),
                JFF = NULL,
                JV = NULL,
                JGG = NULL,
                JW = NULL)
    class(mod) <- "dlm"
    return(mod)
}

###### ARMA
dlmModARMA <- function(ar=NULL, ma=NULL, sigma2=1, dV, m0, C0)
{
    if (is.matrix(sigma2) && (m <- nrow(sigma2)) > 1)
    { ## multivariate
        if (ncol(sigma2) != m)
            stop("sigma2 must be a square matrix")
        r <- max(p <- length(ar), (q <- length(ma)) + 1)
        k <- m * r
        if (hasArg("dV")) {
            if ( length(dV) != m )
                stop("Incompatible dimensions of arguments")
        }
        else
            dV <- rep(0,m)
        if (hasArg("m0")) {
            if ( length(m0) != k )
                stop("Incompatible dimensions of arguments")
        }
        else
            m0 <- rep(0,k)
        if (hasArg("C0")) {
            if ( !( nrow(C0) == k && ncol(C0) == k ))
                stop("Incompatible dimensions of arguments")
        }
        else
            C0 <- 1e7*diag(nrow=k)
        FF <- matrix(0,m,k)
        FF[row(FF) == col(FF)] <- 1
        GG <- matrix(0,k,k)
        for (i in seq(length.out=p))  GG[ (1 + (i-1) * m):(i * m), 1:m ] <- ar[[i]]
        GG[row(GG) == col(GG) - m] <- 1
        R <- matrix(0,k,m)
        R[row(R) == col(R)] <- 1
        for (i in seq(length.out=q)) R[ (1 + i * m):((i+1) * m), ] <- ma[[i]]
        mod <- list(
                    m0 = m0,
                    C0 = C0,
                    FF = FF,
                    V = diag(dV,nrow=m),
                    GG = GG,
                    W = R %*% sigma2 %*% t(R),
                    JFF = NULL,
                    JV = NULL,
                    JGG = NULL,
                    JW = NULL)
    }
    else
    { ## univariate
        r <- max(p <- length(ar), (q <- length(ma)) + 1)
        if (hasArg("dV")) {
            if ( length(dV) != 1 )
                stop("Incompatible dimensions of arguments")
        }
        else
            dV <- 0
        if (hasArg("m0")) {
            if ( length(m0) != r )
                stop("Incompatible dimensions of arguments")
        }
        else
            m0 <- rep(0,r)
        if (hasArg("C0")) {
            if ( !( nrow(C0) == r && ncol(C0) == r ))
                stop("Incompatible dimensions of arguments")
        }
        else
            C0 <- 1e7*diag(nrow=r)
        GG <- matrix(0,r,r)
        if ( p > 0 ) GG[ 1:p, 1 ] <- ar
        GG[row(GG) == col(GG) - 1] <- 1
        R <- rep(0,r)
        R[1] <- 1
        if ( q > 0 ) R[ 1 + 1:q ] <- ma
        mod <- list(
                    m0 = m0,
                    C0 = C0,
                    FF = matrix(c(1,rep(0,r-1)),nrow=1),
                    V = matrix(dV),
                    GG = GG,
                    W = sigma2 * crossprod(matrix(R,nrow=1)),
                    JFF = NULL,
                    JV = NULL,
                    JGG = NULL,
                    JW = NULL)
    }
    class(mod) <- "dlm"
    return(mod)
}


###### Addition for "dlm" objects
"+.dlm" <- function(mod1, mod2)
{
    if ( (m <- nrow(mod1$FF)) != nrow(mod2$FF) )
        stop("Incompatible models")
    if ( (x1 <- !is.null(mod1$X)) | (x2 <- !is.null(mod2$X)) ) {
        if ( x1 && x2 && ( nrow(mod1$X) != nrow(mod2$X) ) )
            stop("Number of rows of mod1$X and mod2$X must match")
        plus <- function(u,v) ifelse( u==0, 0, u+v )
        p1 <- ncol(mod1$FF)
        p2 <- ncol(mod2$FF)
        r <- if ( x1 ) ncol(mod1$X) else 0
        if ( (one <- !is.null(mod1$JFF)) | (two <- !is.null(mod2$JFF)) ) {
            JFF1 <- if ( one ) mod1$JFF else matrix(0,m,p1)
            JFF2 <- if ( two ) plus(mod2$JFF,r) else matrix(0,m,p2)
            JFF <- cbind(JFF1, JFF2)
        } else
        JFF <- NULL
        if ( (one <- !is.null(mod1$JGG)) | (two <- !is.null(mod2$JGG)) ) {
            JGG1 <- if ( one ) mod1$JGG else matrix(0,p1,p1)
            JGG2 <- if ( two ) plus(mod2$JGG,r) else matrix(0,p2,p2)
            JGG <- bdiag(list(JGG1, JGG2))
        } else
        JGG <- NULL
        if ( (one <- !is.null(mod1$JW)) | (two <- !is.null(mod2$JW)) ) {
            JW1 <- if ( one ) mod1$JW else matrix(0,p1,p1)
            JW2 <- if ( two ) plus(mod2$JW,r) else matrix(0,p2,p2)
            JW <- bdiag(list(JW1, JW2))
        } else
        JW <- NULL
        if ( (one <- !is.null(mod1$JV)) | (two <- !is.null(mod2$JV)) ) {
            if ( one && two )
                stop("mod1$V and mod2$V cannot be both time-varying")
            if ( one ) {
                if ( any(mod2$V != 0) ) {
                    mod2$V[] <- 0
                    warning("the value of mod2$V has been discarded")
                }
                JV <- mod1$JV
            }
            if ( two ) {
                if ( any(mod1$V != 0) ) {
                    mod1$V[] <- 0
                    warning("the value of mod1$V has been discarded")
                }
                JV <- plus(mod2$JV,r)
            }
        } else
        JV <- NULL
        mod <- list(
                    m0 = c(mod1$m0, mod2$m0),
                    C0 = bdiag(list(mod1$C0, mod2$C0)),
                    FF = cbind(mod1$FF, mod2$FF),
                    V = mod1$V + mod2$V,
                    GG = bdiag(list(mod1$GG, mod2$GG)),
                    W = bdiag(list(mod1$W, mod2$W)),
                    JFF = JFF,
                    JV = JV,
                    JGG = JGG,
                    JW = JW,
                    X = switch( x1 + 2 * x2,
                    mod1$X,
                    mod2$X,
                    cbind( mod1$X, mod2$X) ))
    } else
    mod <- list(
                m0 = c(mod1$m0, mod2$m0),
                C0 = bdiag(list(mod1$C0, mod2$C0)),
                FF = cbind(mod1$FF, mod2$FF),
                V = mod1$V + mod2$V,
                GG = bdiag(list(mod1$GG, mod2$GG)),
                W = bdiag(list(mod1$W, mod2$W)),
                JFF = NULL,
                JV = NULL,
                JGG = NULL,
                JW = NULL)
    class(mod) <- "dlm"
    return(mod)
}


###### Outer sum of DLMs
"dlmSum" <- function(...)
{
    if ( (narg <- nargs()) == 1  ) {
        if ( is.dlm(...) )
            return(...)
        else
            stop("Argument must be a \"dlm\" object")
    }
    else
    {
        args <- list(...)
        if ( !all(sapply(args, is.dlm)) )
            stop("Arguments must be \"dlm\" objects")
        if ( narg > 2 )
            return( dlmSum(args[[1]], do.call("dlmSum", args[-1])))
        ## now we are in the case narg == 2
        if ((x1 <- !is.null(args[[1]]$X)) | (x2 <- !is.null(args[[2]]$X)))
            stop("Sum of dlm's is only implemented for constant models")
        else
            mod <- list(m0 = c(args[[1]]$m0, args[[2]]$m0),
                        C0 = bdiag(list(args[[1]]$C0, args[[2]]$C0)),
                        FF = bdiag(list(args[[1]]$FF, args[[2]]$FF)),
                        V = bdiag(list(args[[1]]$V, args[[2]]$V)),
                        GG = bdiag(list(args[[1]]$GG, args[[2]]$GG)),
                        W = bdiag(list(args[[1]]$W, args[[2]]$W)),
                        JFF = NULL, JV = NULL, JGG = NULL, JW = NULL)
        class(mod) <- "dlm"
        return(mod)
    }
}

"%+%" <- function(x, y)
    dlmSum(x, y)


###### Print method for "dlm" objects
print.dlm <- function(x, ...) {
    nmRef <- c("FF","V","GG","W","JFF","JV","JGG","JW","X","m0","C0")
    nm <- names(x)
    what <- !sapply(x, is.null)
    ind <- match(nmRef,nm)
    ind <- ind[!is.na(ind)]
    good <- nm[ind][what[ind]]
    if ( !is.na(match("X",good)) && nrow(x$X) > 2 ) {
        x$X <- rbind(formatC(x$X[1:2,,drop=FALSE],digits=4,format="fg"),
                       c("...",rep("",ncol(x$X)-1)))
    }
    print(x[match(good,nm)], quote=FALSE)
    invisible(x)
}

###### Negative loglikelihood
dlmLL <- function(y, mod, debug=FALSE)
{
    ## calculations based on singular value decomposition
    ## Note: V must be nonsingular
    ## The C code relies on the order of the elements in 'mod'
    ## mod = list(m0, C0, FF, V, GG, W)
     if (!debug) {
        storage.mode(y) <- "double"
        matchnames <- match(c("m0", "C0", "FF", "V", "GG", "W"), names(mod))
        for (i in matchnames)
            storage.mode(mod[[i]]) <- "double"
        ## define flags for time-varying components
        if (is.null(mod$JFF))
            tvFF <- FALSE
        else {
            tvFF <- TRUE
            nz <- mod$JFF != 0
            mod$JFF <- cbind(row(mod$JFF)[nz], col(mod$JFF)[nz], mod$JFF[nz]) - 1
            storage.mode(mod$JFF) <- "integer"
        }
        if (is.null(mod$JV))
            tvV <- FALSE
        else {
            tvV <- TRUE
            nz <- mod$JV != 0
            mod$JV <- cbind(row(mod$JV)[nz], col(mod$JV)[nz], mod$JV[nz]) - 1
            storage.mode(mod$JV) <- "integer"
        }
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz]) - 1
            storage.mode(mod$JGG) <- "integer"
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz]) - 1
            storage.mode(mod$JW) <- "integer"
        }
        if (any(c(tvFF,tvV,tvGG,tvW))) {
            mod <- mod[match(c("m0", "C0", "FF", "V", "GG", "W",
                               "JFF", "JV", "JGG", "JW", "X"), names(mod))]
            return(.Call("dlmLL", y, mod, tvFF, tvV, tvGG, tvW, PACKAGE="dlm"))
        } else {
            mod <- mod[match(c("m0", "C0", "FF", "V", "GG", "W"), names(mod))]
            return(.Call("dlmLL0", y, mod, PACKAGE="dlm"))
        }
    }
    else {
        eps <- .Machine$double.eps^.3
        y <- as.matrix(y)
        n <- nrow(y)
        ll <- 0 # negative loglikelihood
        ## define flags for time-varying components
        if (is.null(mod$JFF))
            tvFF <- FALSE
        else {
            tvFF <- TRUE
            nz <- mod$JFF != 0
            mod$JFF <- cbind(row(mod$JFF)[nz], col(mod$JFF)[nz], mod$JFF[nz])
        }
        if (is.null(mod$JV))
            tvV <- FALSE
        else {
            tvV <- TRUE
            nz <- mod$JV != 0
            mod$JV <- cbind(row(mod$JV)[nz], col(mod$JV)[nz], mod$JV[nz])
        }
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz])
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz])
        }
        tvFV <- tvFF || tvV
        ## preliminary calculations, if possible (non time-varying case)
        if ( !tvV ) {
            tmp <- La.svd(mod$V,nu=0)
            Dv <- sqrt(tmp$d)
            if (any(Dv < eps)) {
              Dv <- pmax(Dv, eps)
              warning("a numerically singular 'V' has been slightly perturbed to make it nonsingular")
            }
            Dv.inv <- 1/Dv
            sqrtVinv <- Dv.inv * tmp$vt # t() %*% () = V^{-1}
            sqrtV <- Dv * tmp$vt # t()%*%() = V
            if ( !tvFF )
                tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
        }
        if ( !tvW ) {
            svdW <- La.svd(mod$W,nu=0)
            sqrtW <- sqrt(svdW$d) * svdW$vt # t()%*%() = W
        }
        tmp <- La.svd(mod$C0,nu=0)
        Ux <- t(tmp$vt); Dx <- sqrt(tmp$d)
        for (i in seq(length = n)) {
            ## set time-varying matrices
            if ( tvFF )
                mod$FF[mod$JFF[,-3,drop=FALSE]] <- mod$X[i,mod$JFF[,3]]
            if ( tvV ) {
                mod$V[mod$JV[,-3,drop=FALSE]] <- mod$X[i,mod$JV[,3]]
                tmp <- La.svd(mod$V,nu=0)
                Dv <- sqrt(tmp$d)
                Dv.inv <- 1/Dv; Dv.inv[abs(Dv.inv)==Inf] <- 0
                sqrtVinv <- Dv.inv * tmp$vt
                sqrtV <- sqrt(tmp$d) * tmp$vt # t()%*%() = V
            }
            if ( tvGG )
                mod$GG[mod$JGG[,-3,drop=FALSE]] <- mod$X[i,mod$JGG[,3]]
            if ( tvW ) {
                mod$W[mod$JW[,-3,drop=FALSE]] <- mod$X[i,mod$JW[,3]]
                svdW <- La.svd(mod$W,nu=0)
                sqrtW <- sqrt(svdW$d) * svdW$vt # t()%*%() = W
            }
            if ( tvFV )
                tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)

            if (!any(whereNA <- is.na(y[i, ]))) { ## No missing values

                ## prior
                a <- mod$GG %*% mod$m0
                tmp <- La.svd(rbind( Dx*t(mod$GG%*%Ux), sqrtW ), nu=0)
                Ux.prior <- t(tmp$vt)
                Dx.prior <- tmp$d
                ## one-step forecast
                f <- mod$FF %*% a
                tmp <- La.svd(rbind( Dx.prior*t(mod$FF%*%Ux.prior), sqrtV ), nu=0)
                Uy <- t(tmp$vt)
                Dy <- tmp$d
                ## posterior
                D.inv <- 1/Dx.prior
                D.inv[abs(D.inv)==Inf] <- 0
                tmp <- La.svd(rbind(sqrtVinv%*%mod$FF%*%Ux.prior,
                                    diag(x=D.inv,nrow=length(D.inv))), nu=0)
                Ux <- Ux.prior %*% t(tmp$vt)
                Dx <- 1/tmp$d
                Dx[abs(Dx)==Inf] <- 0
                e <- as.matrix(y[i,]-f)
                mod$m0 <- a + crossprod(Dx*t(Ux)) %*%
                    tF.Vinv %*% e
                ## update scaled negative loglikelihood
                ll <- ll + 2*sum(log(Dy)) + crossprod(crossprod(Uy,e)/Dy)

            } else {
                if (all(whereNA)) { ## All components missing

                    ## prior & posterior
                    mod$m0 <- mod$GG %*% mod$m0
                    tmp <- La.svd(rbind( Dx*t(mod$GG%*%Ux), sqrtW ), nu=0)
                    Ux <- t(tmp$vt)
                    Dx <- tmp$d

                } else { ## Some components missing

                    good <- !whereNA
                    tmp <- La.svd(mod$V[good, good], nu=0)
                    Dv <- sqrt(tmp$d)
                    Dv.inv <- 1/Dv; Dv.inv[abs(Dv.inv)==Inf] <- 0
                    sqrtVinvTMP <- Dv.inv * tmp$vt
                    tF.VinvTMP <- t(mod$FF[good,,drop=FALSE]) %*% crossprod(sqrtVinvTMP)
                    sqrtVTMP <- Dv * tmp$vt
                    ## prior
                    a <- mod$GG %*% mod$m0
                    tmp <- La.svd(rbind( Dx*t(mod$GG%*%Ux), sqrtW ), nu=0)
                    Ux.prior <- t(tmp$vt)
                    Dx.prior <- tmp$d
                    ## one-step forecast
                    f <- mod$FF[good,,drop=FALSE] %*% a
                    tmp <- La.svd(rbind( Dx.prior*t(mod$FF[good,,drop=FALSE]%*%Ux.prior), sqrtVTMP ), nu=0)
                    Uy <- t(tmp$vt)
                    Dy <- tmp$d
                    ## posterior
                    D.inv <- 1/Dx.prior
                    D.inv[abs(D.inv)==Inf] <- 0
                    tmp <- La.svd(rbind(sqrtVinvTMP%*%mod$FF[good,,drop=FALSE]%*%Ux.prior,
                                        diag(x=D.inv,nrow=length(D.inv))), nu=0)
                    Ux <- Ux.prior %*% t(tmp$vt)
                    Dx <- 1/tmp$d
                    Dx[abs(Dx)==Inf] <- 0
                    e <- as.matrix(y[i,good]-f)
                    mod$m0 <- a + crossprod(Dx*t(Ux)) %*%
                        tF.VinvTMP %*% e
                    ## update scaled negative loglikelihood
                    ll <- ll + 2*sum(log(Dy)) + crossprod(crossprod(Uy,e)/Dy)
                }
            }
        }
    }
    return(drop(ll)*0.5) # negative loglikelihood
}


"dlmMLE" <- function(y, parm, build, method = "L-BFGS-B",
                     ..., debug = FALSE)
{
    logLik <- function(parm, ...)
    {
        mod <- build(parm, ...)
        return(dlmLL(y=y, mod=mod, debug=debug))
    }
    out <- optim(parm, logLik, method=method, ...)
    return(out)
}


dlmFilter <- function(y, mod, debug = FALSE, simplify = FALSE)
{
    ## Note: V must be nonsingular. It will be forced to be so,
    ## with a warning, otherwise (time-invariant case only).
    eps <- .Machine$double.eps^.3
    mod1 <- mod
    yAttr <- attributes(y)
    ytsp <- tsp(y)
    y <- as.matrix(y)
    timeNames <- dimnames(y)[[1]]
    stateNames <- names(mod$m0)
    if (!debug) {
        storage.mode(y) <- "double"
        matchnames <- match(c("m0", "C0", "FF", "V", "GG", "W"), names(mod))
        for (i in matchnames)
            storage.mode(mod[[i]]) <- "double"
        if ("X" %in% names(mod))
          storage.mode(mod$X) <- "double"
        ## define flags for time-varying components
        if (is.null(mod$JFF))
            tvFF <- FALSE
        else {
            tvFF <- TRUE
            nz <- mod$JFF != 0
            mod$JFF <- cbind(row(mod$JFF)[nz], col(mod$JFF)[nz], mod$JFF[nz]) - 1
            storage.mode(mod$JFF) <- "integer"
        }
        if (is.null(mod$JV))
            tvV <- FALSE
        else {
            tvV <- TRUE
            nz <- mod$JV != 0
            mod$JV <- cbind(row(mod$JV)[nz], col(mod$JV)[nz], mod$JV[nz]) - 1
            storage.mode(mod$JV) <- "integer"
        }
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz]) - 1
            storage.mode(mod$JGG) <- "integer"
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz]) - 1
            storage.mode(mod$JW) <- "integer"
        }
        if (any(c(tvFF,tvV,tvGG,tvW))) {
            mod <- mod[match(c("m0", "C0", "FF", "V", "GG", "W",
                               "JFF", "JV", "JGG", "JW", "X"), names(mod))]
            ans <- .Call("dlmFilter", y, mod, tvFF, tvV, tvGG, tvW, PACKAGE = "dlm")
        } else {
            mod <- mod[match(c("m0", "C0", "FF", "V", "GG", "W"), names(mod))]
            ans <- .Call("dlmFilter0", y, mod, PACKAGE = "dlm")
        }
        names(ans) <- c("m", "U.C", "D.C", "a", "U.R", "D.R", "f")
    }
    else {
        ## define flags for time-varying components
        if (is.null(mod$JFF))
            tvFF <- FALSE
        else {
            tvFF <- TRUE
            nz <- mod$JFF != 0
            mod$JFF <- cbind(row(mod$JFF)[nz], col(mod$JFF)[nz], mod$JFF[nz])
        }
        if (is.null(mod$JV))
            tvV <- FALSE
        else {
            tvV <- TRUE
            nz <- mod$JV != 0
            mod$JV <- cbind(row(mod$JV)[nz], col(mod$JV)[nz], mod$JV[nz])
        }
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz])
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz])
        }
        tvFV <- tvFF || tvV
        m <- rbind(mod$m0,matrix(0,nrow=nrow(y),ncol=length(mod$m0))) # filtered values
        a <- matrix(0,nrow=nrow(y),ncol=length(mod$m0))
        f <- matrix(0,nrow=nrow(y),ncol=ncol(y))
        U.C <- vector(1+nrow(y),mode="list")
        D.C <- matrix(0,1+nrow(y),length(mod$m0))
        U.R <- vector(nrow(y),mode="list")
        D.R <- matrix(0,nrow(y),length(mod$m0))
        ## preliminary calculations, if possible (non time-varying case)
        if ( !tvV ) {
            tmp <- La.svd(mod$V,nu=0)
            Uv <- t(tmp$vt); Dv <- sqrt(tmp$d)
            if (any(Dv < eps)) {
              Dv <- pmax(Dv, eps)
              warning("a numerically singular 'V' has been slightly perturbed to make it nonsingular")
            }
            Dv.inv <- 1/Dv
            sqrtVinv <- Dv.inv * tmp$vt # t() %*% () = V^{-1}
            if ( !tvFF )
                tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
        }
        if ( !tvW ) {
            svdW <- La.svd(mod$W,nu=0)
            sqrtW <- sqrt(svdW$d) * svdW$vt # t()%*%() = W
        }
        tmp <- La.svd(mod$C0,nu=0)
        U.C[[1]] <- t(tmp$vt)
        D.C[1,] <- sqrt(tmp$d)
        for (i in seq(length=nrow(y))) {
            ## set time-varying matrices
            if ( tvFF )
                mod$FF[mod$JFF[,-3,drop=FALSE]] <- mod$X[i,mod$JFF[,3]]
            if ( tvV ) {
                mod$V[mod$JV[,-3,drop=FALSE]] <- mod$X[i,mod$JV[,3]]
                tmp <- La.svd(mod$V,nu=0)
                Uv <- t(tmp$vt); Dv <- sqrt(tmp$d)
                Dv.inv <- 1/Dv; Dv.inv[abs(Dv.inv)==Inf] <- 0
                sqrtVinv <- Dv.inv * tmp$vt
            }
            if ( tvGG )
                mod$GG[mod$JGG[,-3,drop=FALSE]] <- mod$X[i,mod$JGG[,3]]
            if ( tvW ) {
                mod$W[mod$JW[,-3,drop=FALSE]] <- mod$X[i,mod$JW[,3]]
                svdW <- La.svd(mod$W,nu=0)
                sqrtW <- sqrt(svdW$d) * svdW$vt # t()%*%() = W
            }
            if ( tvFV )
                tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)

            if (!any(whereNA <- is.na(y[i, ]))) { ## No missing values

                ## prior
                a[i,] <- mod$GG %*% m[i,]
                tmp <- La.svd(rbind( D.C[i,]*t(mod$GG%*%U.C[[i]]), sqrtW ), nu=0)
                U.R[[i]] <- t(tmp$vt)
                D.R[i,] <- tmp$d
                ## one-step forecast
                f[i,] <- mod$FF %*% a[i,]
                ## posterior
                D.Rinv <- 1/D.R[i,]
                D.Rinv[abs(D.Rinv)==Inf] <- 0
                tmp <- La.svd(rbind(sqrtVinv %*% mod$FF %*% U.R[[i]],
                                    diag(x=D.Rinv,nrow=length(D.Rinv))), nu=0)
                U.C[[i+1]] <- U.R[[i]] %*% t(tmp$vt)
                foo <- 1/tmp$d; foo[abs(foo)==Inf] <- 0
                D.C[i+1,] <- foo
                m[i+1,] <- a[i,] + crossprod(D.C[i+1,]*t(U.C[[i+1]])) %*%
                    tF.Vinv %*% as.matrix(y[i,]-f[i,])

            } else {
                if (all(whereNA)) { ## All components missing

                    ## prior & posterior
                    m[i+1,] <- a[i,] <- mod$GG %*% m[i,]
                    tmp <- La.svd(rbind( D.C[i,]*t(mod$GG%*%U.C[[i]]), sqrtW ), nu=0)
                    U.C[[i+1]] <- U.R[[i]] <- t(tmp$vt)
                    D.C[i+1,] <- D.R[i,] <- tmp$d
                    ## one-step forecast
                    f[i,] <- mod$FF %*% a[i,]

                } else { ## Some components missing

                    good <- !whereNA
                    tmp <- La.svd(mod$V[good, good], nu=0)
                    Dv <- sqrt(tmp$d)
                    Dv.inv <- 1/Dv; Dv.inv[abs(Dv.inv)==Inf] <- 0
                    sqrtVinvTMP <- Dv.inv * tmp$vt
                    tF.VinvTMP <- t(mod$FF[good,,drop=FALSE]) %*% crossprod(sqrtVinvTMP)
                    ## prior
                    a[i,] <- mod$GG %*% m[i,]
                    tmp <- La.svd(rbind( D.C[i,]*t(mod$GG%*%U.C[[i]]), sqrtW ), nu=0)
                    U.R[[i]] <- t(tmp$vt)
                    D.R[i,] <- tmp$d
                    ## one-step forecast
                    f[i,] <- mod$FF %*% a[i,]
                    ## posterior
                    D.Rinv <- 1/D.R[i,]
                    D.Rinv[abs(D.Rinv)==Inf] <- 0
                    tmp <- La.svd(rbind(sqrtVinvTMP %*% mod$FF[good,,drop=FALSE] %*% U.R[[i]],
                                        diag(x=D.Rinv,nrow=length(D.Rinv))), nu=0)
                    U.C[[i+1]] <- U.R[[i]] %*% t(tmp$vt)
                    foo <- 1/tmp$d; foo[abs(foo)==Inf] <- 0
                    D.C[i+1,] <- foo
                    m[i+1,] <- a[i,] + crossprod(D.C[i+1,]*t(U.C[[i+1]])) %*%
                        tF.VinvTMP %*% as.matrix(y[i,good]-f[i,good])
                }
            }
        }

        ans <- list(m=m,U.C=U.C,D.C=D.C,a=a,U.R=U.R,D.R=D.R,f=f)
    }

    ans$m <- drop(ans$m); ans$a <- drop(ans$a); ans$f <- drop(ans$f)
    attributes(ans$f) <- yAttr
    if (!is.null(ytsp)) {
        tsp(ans$a) <- ytsp
        tsp(ans$m) <- c(ytsp[1] - 1/ytsp[3], ytsp[2:3])
        class(ans$a) <- class(ans$m) <- if (length(mod$m0) > 1) c("mts","ts") else "ts"
    }
    if (!(is.null(timeNames) && is.null(stateNames)))
        if (is.matrix(ans$a))
        {
            dimnames(ans$a) <- list(timeNames, stateNames)
            dimnames(ans$m) <- list(if(is.null(timeNames)) NULL else c("",timeNames),
                                    stateNames)
        }
        else
            if (!is.null(timeNames))
            {
                names(ans$a) <- timeNames
                names(ans$m) <- c("", timeNames)
            }
    if (simplify)
        ans <- c(mod=list(mod1), ans)
    else
    {
        attributes(y) <- yAttr
        ans <- c(y=list(y), mod=list(mod1), ans)
    }
    class(ans) <- "dlmFiltered"
    return(ans)
}

dlmSmooth <- function(y, ...)
    UseMethod("dlmSmooth", y)

dlmSmooth.default <- function(y, mod, ...)
{
    dlmSmooth(dlmFilter(y, mod, ...), ...)
}

dlmSmooth.dlmFiltered <- function(y, ..., debug = FALSE)
{
    big <- 1 / sqrt(.Machine$double.eps)
    mod <- c(y[match(c("m", "U.C", "D.C", "a", "U.R", "D.R"),names(y))],
             y$mod[match(c("GG", "W", "JGG", "JW", "X"), names(y$mod))])
    mAttr <- attributes(mod$m)
    mod$m <- as.matrix(mod$m)
    mod$a <- as.matrix(mod$a)
    if (!debug) {
        ## define flags for time-varying components
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz]) - 1
            storage.mode(mod$JGG) <- "integer"
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz]) - 1
            storage.mode(mod$JW) <- "integer"
        }
        if (tvGG || tvW) {
            ans <- .Call("dlmSmooth", mod, tvGG, tvW, big, PACKAGE="dlm")
            names(ans) <- c("s", "U.S", "D.S")
        }
        else {
            ans <- .Call("dlmSmooth0", mod, big, PACKAGE="dlm")
            names(ans) <- c("s", "U.S", "D.S")
        }
        } else {
        if (is.null(mod$JGG))
            tvGG <- FALSE
        else {
            tvGG <- TRUE
            nz <- mod$JGG != 0
            mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz])
        }
        if (is.null(mod$JW))
            tvW <- FALSE
        else {
            tvW <- TRUE
            nz <- mod$JW != 0
            mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz])
        }
        n <- length(mod$U.R) # number of obs
        p <- NCOL(mod$m) # dimension of state vector
        s <- rbind(matrix(0,n,p), mod$m[n+1,])
        U.S <- vector("list", length=n+1)
        U.S[[n+1]] <- mod$U.C[[n+1]]
        D.S <- rbind(matrix(0,n,p), mod$D.C[n+1,])
        ## preliminary calculations, if possible (time-invariant case)
        if ( !tvW ) {
            tmp <- La.svd(mod$W,nu=0)
            Dw <- sqrt(tmp$d)
            Dw.inv <- pmin(1/Dw, big)
            sqrtWinv <- Dw.inv * tmp$vt # t()%*%() = W^(-1)
        }
        if (n > 0)
            for (i in n:1)
            {
                ## set relevant time-varying matrices
                if ( tvGG )
                    mod$GG[mod$JGG[,-3,drop=FALSE]] <- mod$X[i,mod$JGG[,3]]
                if ( tvW ) {
                    mod$W[mod$JW[,-3,drop=FALSE]] <- mod$X[i,mod$JW[,3]]
                    tmp <- La.svd(mod$W,nu=0)
                    Dw <- sqrt(tmp$d)
                    Dw.inv <- pmin(1/Dw, big)
                    sqrtWinv <- Dw.inv * tmp$vt # t()%*%() = W^(-1)
                }
                Dinv <- 1/mod$D.R[i,]; Dinv[abs(Dinv)==Inf] <- 0
                H <- crossprod(mod$D.C[i,]*t(mod$U.C[[i]])) %*%
                    t(mod$GG) %*% crossprod(Dinv*t(mod$U.R[[i]]))
                Dinv <- 1/mod$D.C[i,]
                Dinv[abs(Dinv)==Inf] <- 0
                tmp <- La.svd(rbind( sqrtWinv%*%mod$GG, Dinv*t(mod$U.C[[i]])), nu=0)
                Dinv <- 1/tmp$d
                Dinv[abs(Dinv)==Inf] <- 0
                tmp <- La.svd(rbind(Dinv*tmp$vt, D.S[i+1,]*t(H%*%U.S[[i+1]])))
                U.S[[i]] <- t(tmp$vt)
                D.S[i,] <- tmp$d
                s[i,] <- mod$m[i,] + H %*% (s[i+1,]-mod$a[i,])
            }
        ans <- list(s=s, U.S=U.S, D.S=D.S)
    }
    attributes(ans$s) <- mAttr
                                        #     if (!is.null(tsp(mod$m))) {
#         tsp(ans$s) <- tsp(mod$m)
#         class(ans$s) <- if (NCOL(ans$s) > 1) c("mts","ts") else "ts"
#     }
#     dimnames(ans$s) <- dimnames(mod$m)
    return(ans)
}


dlmBSample <- function(modFilt)
{
    eps <- .Machine$double.eps^.4
    mod <- c(modFilt[match(c("m", "U.C", "D.C", "a", "U.R", "D.R"),names(modFilt))],
             modFilt$mod[match(c("GG", "W", "JGG", "JW", "X"), names(modFilt$mod))])
    n <- length(mod$U.R) # number of obs
    p <- NCOL(mod$m) # dimension of state vector
    mtsp <- tsp(mod$m)
    if (p==1) {
        dim(mod$m) <- c(n+1, 1)
        dim(mod$a) <- c(n, 1)
    }
    if (is.null(mod$JGG))
        tvGG <- FALSE
    else {
        tvGG <- TRUE
        nz <- mod$JGG != 0
        mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz])
    }
    if (is.null(mod$JW))
        tvW <- FALSE
    else {
        tvW <- TRUE
        nz <- mod$JW != 0
        mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz])
    }
    tvGW <- tvGG || tvW
    theta <- matrix(0,n+1,p)
    ## preliminary calculations, if possible (time-invariant case)
    if ( !tvW ) {
        tmp <- La.svd(mod$W,nu=0)
        Dw <- sqrt(tmp$d)
        Dw <- pmax(Dw, eps)
        Dw.inv <- 1/Dw
        sqrtWinv <- Dw.inv * tmp$vt # t()%*%() = W^(-1)
        if ( !tvGG ) tG.Winv <- t(mod$GG) %*% crossprod(sqrtWinv)
    }
    ## generate last theta
    theta[n+1,] <- mod$m[n+1,] + mod$U.C[[n+1]] %*% matrix(mod$D.C[n+1,]*rnorm(p))
    ## generate all the other theta's
    for (i in (n:1))
    {
        if ( tvGG )
            mod$GG[mod$JGG[,-3,drop=FALSE]] <- mod$X[i,mod$JGG[,3]]
        if ( tvW ) {
            mod$W[mod$JW[,-3,drop=FALSE]] <- mod$X[i,mod$JW[,3]]
            tmp <- La.svd(mod$W,nu=0)
            Dw <- sqrt(tmp$d)
            Dw <- pmax(Dw, eps)
            Dw.inv <- 1/Dw
            sqrtWinv <- Dw.inv * tmp$vt # t()%*%() = W^(-1)
        }
        if ( tvGW )
            tG.Winv <- t(mod$GG) %*% crossprod(sqrtWinv)
        D.inv <- 1/mod$D.C[i,]; D.inv[abs(D.inv)==Inf] <- 0
        tmp <- La.svd(rbind(sqrtWinv %*% mod$GG %*% mod$U.C[[i]],
                            diag(x=D.inv,nrow=length(D.inv))), nu=0)
        U.H <- mod$U.C[[i]] %*% t(tmp$vt)
        D.H <- 1/tmp$d; D.H[abs(D.H)==Inf] <- 0
        h <- mod$m[i,] + crossprod(D.H*t(U.H)) %*%
            tG.Winv %*% (t(theta[i+1,,drop=F])-mod$a[i,])
        theta[i,] <- h + U.H %*% matrix(D.H*rnorm(p))
    }
    if (!is.null(mtsp))
    {
        theta <- drop(theta)
        tsp(theta) <- mtsp
        class(theta) <- if (p > 1) c("mts","ts") else "ts"
    }
    return(theta=theta)
}


rwishart <- function(df, p = nrow(SqrtSigma), Sigma,
                     SqrtSigma = diag(p))
{
    ## generate a Wishart-distributed matrix - from S-news, due to B. Venables
    ## note: Sigma = crossprod(SqrtSigma), and df must be integer
    if (!missing(Sigma))
    {
        ## compute SqrtSigma
        tmp <- svd(Sigma)
        SqrtSigma <- sqrt(tmp$d) * t(tmp$u)
    }
    if((Ident <- missing(SqrtSigma)) && missing(p))
        stop("either p, Sigma or SqrtSigma must be specified")
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, df:(df-p+1)))
    if(p > 1)
    {
        pseq <- 1:(p-1)
        Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
    }
    if(Ident)
        crossprod(Z)
    else
        crossprod(Z %*% SqrtSigma)
}


dlmSvd2var <- function(u, d)
{
    if (is.matrix(u)) {
        if ( nrow(u) != length(d) || ncol(u) != (n <- length(d)) )
            stop("inconsistent dimensions")
        return(tcrossprod(rep(d, each = n) * u))
    }
    if (is.list(u)) {
        if ( length(u) != NROW(d) )
            stop("length of 'u' must be equal to the number of rows of 'd'")
        n <- NCOL(d)
        return(lapply(seq(along = u), function(i)
                      tcrossprod(rep(d[i,], each = n) * u[[i]])))

    }
    stop("wrong argument 'u'")
}

dlmForecast <- function(mod, nAhead=1, method=c("plain","svd"), sampleNew=FALSE) {
    method <- match.arg(method)

    if (is.dlmFiltered(mod)) {
        modFuture <- mod$mod
        lastObsIndex <- NROW(mod$m)
        modFuture$C0 <- with(mod, dlmSvd2var(U.C[[lastObsIndex]], D.C[lastObsIndex,]))
        if (is.ts(mod$m))
            modFuture$m0 <- window(mod$m, start=end(mod$m))
        else {
            modFuture$m0 <- window(mod$m, start=lastObsIndex)
            tsp(modFuture$m0) <- NULL
        }
        mod <- modFuture
    }

    if (! (is.null(mod$JFF) && is.null(mod$JV) && is.null(mod$JGG) && is.null(mod$JW)))
        stop("dlmForecast only works with constant models")

    ytsp <- tsp(mod$m0)
    p <- length(mod$m0)
    m <- nrow(mod$FF)
    a <- rbind(mod$m0, matrix(0,nAhead,p))
    R <- vector("list",nAhead+1)
    R[[1]] <- mod$C0
    f <- matrix(0,nAhead,m)
    Q <- vector("list", nAhead)
    for (it in 1:nAhead) {
        a[it+1,] <- mod$GG %*% a[it,]
        R[[it+1]] <- mod$GG %*% R[[it]] %*% t(mod$GG) + mod$W
        f[it,] <- mod$FF %*% a[it+1,]
        Q[[it]] <- mod$FF %*% R[[it+1]] %*% t(mod$FF) + mod$V
    }
    a <- a[-1,,drop=FALSE]
    R <- R[-1]
    if ( sampleNew ) {
        newStates <- vector("list", sampleNew)
        newObs <- vector("list", sampleNew)
        newS <- matrix(0, nAhead, p)
        newO <- matrix(0, nAhead, m)
        tmp <- La.svd(mod$V,nu=0)
        Ut.V <- tmp$vt; D.V <- sqrt(tmp$d)
        tmp <- La.svd(mod$W,nu=0)
        Ut.W <- tmp$vt; D.W <- sqrt(tmp$d)
        for (i in 1:sampleNew) {
            tmp <- La.svd(R[[1]],nu=0)
            newS[1,] <- crossprod(tmp$vt, rnorm(p, sd=sqrt(tmp$d))) + a[1,]
            newO[1,] <- crossprod(Ut.V, rnorm(m, sd=D.V)) + mod$FF %*% newS[1,]
            if ( nAhead > 1 )
                for (it in 2:nAhead) {
                    newS[it,] <- crossprod(Ut.W, rnorm(p, sd=D.W)) + mod$GG %*% newS[it-1,]
                    newO[it,] <- crossprod(Ut.V, rnorm(m, sd=D.V)) + mod$FF %*% newS[it,]
                }
            newStates[[i]] <- newS
            newObs[[i]] <- newO
        }
        if (!is.null(ytsp)) {
            a <- ts(a, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3])
            f <- ts(f, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3])
            newStates <- lapply(newStates, function(x)
                                ts(x, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3]))
            newObs <- lapply(newObs, function(x)
                             ts(x, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3]))
        }
        ans <- list(a=a, R=R, f=f, Q=Q, newStates=newStates, newObs=newObs)
    } else {
        if (!is.null(ytsp)) {
            a <- ts(a, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3])
            f <- ts(f, start = ytsp[2] + 1/ytsp[3], frequency = ytsp[3])
        }
        ans <- list(a=a, R=R, f=f, Q=Q)
    }


    return(ans)
}

###### One-step forecast errors
residuals.dlmFiltered <- function(object, ...,
                                  type=c("standardized", "raw"), sd=TRUE) {
    if (is.null(object$y))
        stop("\'object\' argument has no \'y\' component")
    type <- match.arg(type)
    if (is.null(object$mod$JFF)) tvFF <- FALSE else tvFF <- TRUE
    if (is.null(object$mod$JV)) tvV <- FALSE else tvV <- TRUE
    FF <- object$mod$FF
    if (!( tvFF || tvV )) { ## constant model
        f <- object$a %*% t(FF)
        res <- drop(object$y - f) # one-step forecasting errors
        if (sd || (type == "standardized")) {
            V <- object$mod$V
            SD <- drop(t(sqrt(sapply(seq(along=object$U.R),
                                     function(i)
                                     diag(crossprod(object$D.R[i,] *
                                                    t(FF%*%object$U.R[[i]])) + V)))))
        }
    } else
    if ( !tvFF ) { ## only V time-varying
        f <- object$a %*% t(FF)
        res <- drop(object$y - f) # one-step forecasting errors
        if (sd || (type == "standardized")) {
            nz <- object$mod$JV != 0
            JV <- cbind(row(object$mod$JV)[nz], col(object$mod$JV)[nz],
                        object$mod$JV[nz])
            V <- object$mod$V
            getSD <- function(i) {
                V[JV[,-3,drop=FALSE]] <- object$mod$X[i,JV[,3]]
                diag(crossprod(object$D.R[i,] * t(FF%*%object$U.R[[i]])) + V)
            }
            SD <- drop(t(sqrt(sapply(seq(along=object$U.R), getSD))))
        }
    } else
    if ( !tvV ) { ## only FF time-varying
        if (!(sd || (type == "standardized"))) {
            nz <- object$mod$JFF != 0
            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],
                         object$mod$JFF[nz])
            getFore <- function(i) {
                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]
                FF %*% object$a[i,]
            }
            f <- drop(t(sapply(seq(along=object$U.R), getFore)))
            res <- drop(object$y - f) # one-step forecasting errors
        } else {
            nz <- object$mod$JFF != 0
            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],
                         object$mod$JFF[nz])
            V <- object$mod$V
            getBoth <- function(i) {
                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]
                c(FF %*% object$a[i,],
                  diag(crossprod(object$D.R[i,] * t(FF%*%object$U.R[[i]])) + V))
            }
            tmp <- t(sapply(seq(along=object$U.R), getBoth))
            m <- ncol(tmp) / 2
            res <- drop(object$y - tmp[,1:m])
            SD <- drop(sqrt(tmp[,-(1:m)]))
        }
    } else { ## both FF and V time-varying
        if (!(sd || (type == "standardized"))) {
            nz <- object$mod$JFF != 0
            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],
                         object$mod$JFF[nz])
            getFore <- function(i) {
                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]
                FF %*% object$a[i,]
            }
            f <- drop(t(sapply(seq(along=object$U.R), getFore)))
            res <- drop(object$y - f) # one-step forecasting errors
        } else {
            nz <- object$mod$JFF != 0
            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],
                         object$mod$JFF[nz])
            nz <- object$mod$JV != 0
            JV <- cbind(row(object$mod$JV)[nz], col(object$mod$JV)[nz],
                        object$mod$JV[nz])
            V <- object$mod$V
            getBoth <- function(i) {
                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]
                V[JV[,-3,drop=FALSE]] <- object$mod$X[i,JV[,3]]
                c(FF %*% object$a[i,],
                  diag(crossprod(object$D.R[i,] * t(FF%*%object$U.R[[i]])) + V))
            }
            tmp <- t(sapply(seq(along=object$U.R), getBoth))
            m <- ncol(tmp) / 2
            res <- drop(object$y - tmp[,1:m])
            SD <- drop(sqrt(tmp[,-(1:m)]))
        }
    }

    if ( type == "standardized" )
        res <- res / SD
    if (sd) {
        if (is.ts(res)) attributes(SD) <- attributes(res) # makes a time series of SD
        return(list(res=res, sd=SD))
    } else
    return(res)
}

###### Diagnostic plots
tsdiag.dlmFiltered <- function (object, gof.lag = 10, ...) {
    stdres <- residuals(object, sd=FALSE)
    if ((ns <- NCOL(stdres)) == 1) {
        oldpar <- par(mfrow = c(3, 1))
        on.exit(par(oldpar))
        plot(stdres, type = "h", main = "Standardized Residuals", ylab = "")
        abline(h = 0)
        acf(stdres, plot = TRUE, main = "ACF of Residuals",
            na.action = na.pass)
        nlag <- gof.lag
        pval <- numeric(nlag)
        for (i in 1:nlag) pval[i] <- Box.test(stdres, i, type = "Ljung-Box")$p.value
        plot(1:nlag, pval, xlab = "lag", ylab = "p value",
             ylim = c(0, 1), main = "p values for Ljung-Box statistic")
        abline(h = 0.05, lty = 2, col = "blue")
    } else {
        ask <- dev.interactive()
        oldpar <- par(mfrow = c(3, 1), oma=c(0, 0, 2, 0), "ask")
        on.exit(par(oldpar))
        hasNames <- !is.null(nm <- attr(stdres,"dimnames")[[2]])
        for (j in 1:ns) {
            plot(stdres[,j], type = "h", main = "Standardized Residuals", ylab = "")
            abline(h = 0)
            if (ask) par(ask=FALSE)
            acf(stdres[,j], plot = TRUE, main = "ACF of Residuals",
                na.action = na.pass)
            nlag <- gof.lag
            pval <- numeric(nlag)
            for (i in 1:nlag) pval[i] <- Box.test(stdres[,j], i, type = "Ljung-Box")$p.value
            plot(1:nlag, pval, xlab = "lag", ylab = "p value",
                 ylim = c(0, 1), main = "p values for Ljung-Box statistic")
            abline(h = 0.05, lty = 2, col = "blue")
            mtext(if (hasNames) nm[j] else paste("Series",j), line=1, outer=TRUE)
            if (ask) par(ask=TRUE)
        }
    }
}

###### Generating a random DLM
dlmRandom <- function(m, p, nobs = 0, JFF, JV, JGG, JW)
{
    ### Assume for now that each of FF, V, GG, W is either fixed or
    ### time-varying in every entry
    FF <- matrix(rnorm(m*p),m,p)
    V <- rwishart(2*m, m)
    GGtrial <- matrix(rnorm(p*p),p,p)
    e <- eigen(GGtrial)
    if ((ab <- max(abs(e$values))) > 1) {
        r <- runif(1)
        GG <- with(e, Re(vectors %*% (r * values / ab * solve(vectors))))
    }
    else
        GG <- GGtrial
    W <- rwishart(2*p, p)
    m0 <- rep(0,p)
    C0 <- diag(nrow = p) * 100
    if (nobs > 0) {
        count <- 0
        if (hasArg(JGG) && JGG == TRUE) {
            tvGG <- TRUE
            JGG <- structure(1:(p*p), dim=dim(GG))
            count <- p * p
        }
        else {
            tvGG <- FALSE
            JGG <- NULL
        }
        if (hasArg(JW) && JW == TRUE) {
            tvW <- TRUE
            JW <- structure(count + 1:(p*p), dim=dim(W))
            count <- count + p * p
        }
        else {
            tvW <- FALSE
            JW <- NULL
            tmp <- La.svd(W,nu=0)
            Ut.W <- tmp$vt; D.W <- sqrt(tmp$d)
        }
        if (hasArg(JFF) && JFF == TRUE) {
            tvFF <- TRUE
            JFF <- structure(count + 1:(m*p), dim=dim(FF))
            count <- count + m * p
        }
        else {
            tvFF <- FALSE
            JFF <- NULL
        }
        if (hasArg(JV) && JV == TRUE) {
            tvV <- TRUE
            JV <- structure(count + 1:(m*m), dim=dim(V))
        }
        else {
            tvV <- FALSE
            JV <- NULL
            tmp <- La.svd(V,nu=0)
            Ut.V <- tmp$vt; D.V <- sqrt(tmp$d)
        }
        if (any(c(tvFF, tvV, tvGG, tvW))) {
            newStates <- matrix(0, nrow = nobs + 1, ncol = p)
            newObs <- matrix(0, nrow = nobs, ncol = m)
            X <- matrix(0, nrow = nobs,
                        ncol = m * p * tvFF + m * m * tvV +
                        p * p * tvGG + p * p * tvW)
            tmp <- La.svd(C0,nu=0)
            newStates[1,] <- crossprod(tmp$vt, rnorm(p, sd=sqrt(tmp$d)))
            for (it in 1:nobs) {
                count <- 0
                if (tvGG) {
                    GGtrial <- matrix(rnorm(p*p),p,p)
                    e <- eigen(GGtrial)
                    if ((ab <- max(abs(e$values))) > 1) {
                        r <- runif(1)
                        GG <- with(e, Re(vectors %*% (r * values / ab * solve(vectors))))
                    }
                    else
                        GG <- GGtrial
                    X[it, 1:(p*p)] <- GG
                    count <- p * p
                }
                if (tvW) {
                    W <- rwishart(2*p, p)
                    X[it, count + 1:(p*p)] <- W
                    count <- count + p * p
                    tmp <- La.svd(W,nu=0)
                    Ut.W <- tmp$vt; D.W <- sqrt(tmp$d)
                }
                if (tvFF) {
                    FF <- matrix(rnorm(m*p),m,p)
                    X[it, count + 1:(m*p)] <- FF
                    count <- count + m * p
                }
                if (tvV) {
                    V <- rwishart(2*m, m)
                    X[it, count + 1:(m*m)] <- V
                    tmp <- La.svd(V,nu=0)
                    Ut.V <- tmp$vt; D.V <- sqrt(tmp$d)
                }
                newStates[it + 1,] <- crossprod(Ut.W, rnorm(p, sd=D.W)) + GG %*% newStates[it,]
                newObs[it,] <- crossprod(Ut.V, rnorm(m, sd=D.V)) + FF %*% newStates[it+1,]
            }
            mod <- dlm(list(m0=m0, C0=C0, FF=FF, V=V, GG=GG, W=W, JFF=JFF, JV=JV,
                            JGG=JGG, JW=JW, X=X))
            ans <- list(mod = mod, theta = newStates[-1,,drop=FALSE], y = newObs)
        }
        else {
            mod <- dlm(list(m0=m0, C0=C0, FF=FF, V=V, GG=GG, W=W))
            tmp <- dlmForecast(mod, nAhead = nobs, sampleNew = 1)
            ans <- list(mod = mod, theta = tmp$newStates[[1]], y = tmp$newObs[[1]])
        }
    }
    else
        ans <- dlm(list(m0=m0, C0=C0, FF=FF, V=V, GG=GG, W=W))
    return(ans)
}

###### Utility functions for MCMC output analysis
mcmcSD <- function(x) {
    ## Monte Carlo standard deviation using Sokal's method
    msg <- "Numeric vector or matrix argument needed"
    if (!is.numeric(x))
        stop(msg)
    univariate <- function(x) {
        l <- floor(30*log10(length(x)))
        ac <- drop(acf(x, lag.max=l, plot=FALSE)$acf)
        tau <- cumsum(c(1, 2 * ac[-1]))
        k <- 0
        while ( (k < 3 * tau[k+1]) &&  (k < l)) k <- k+1
        tau <- tau[k+1]
        if (tau < 0)
            {
                tau <- 1
                warning("estimated s.d. may not be reliable")
            }
        return(sqrt(var(x) * tau / length(x)))
    }
    if (is.null(dim(x)))
        ans <- univariate(x)
    else
        if (is.matrix(x))
            ans <- apply(x, 2, univariate)
        else
            stop(msg)
    return(ans)
}

mcmcMean <- function(x, sd = TRUE)
{
    ## Ergodic means of mcmc output
    ## Estimated standard deviations of means using Sokal's method
    msg <- "Numeric vector or matrix argument needed"
    if (!is.numeric(x))
        stop(msg)
    if (v <- is.null(dim(x))) # univariate input
    {
        mn <- mean(x)
        nm <- deparse(substitute(x))
    }
    else
        if (is.matrix(x)) # multivariate input
        {
            mn <- colMeans(x)
            nm <- colnames(x)
            if ( is.null(nm) )
                nm <- paste(deparse(substitute(x)), 1:NCOL(x), sep='.')
        }
        else
            stop(msg)
    if (sd) {
        univariateSD <- function(x) {
            l <- floor(30*log10(length(x)))
            ac <- acf(x, lag.max=l, plot=FALSE)$acf
            tau <- cumsum(c(1, 2 * ac[-1]))
            k <- 0
            while ( (k < 3 * tau[k+1]) &&  (k < l)) k <- k+1
            tau <- tau[k+1]
            return(sqrt(var(x) * tau / length(x)))
        }
        if ( v )
            sd <- univariateSD(x)
        else
            sd <- apply(x, 2, univariateSD)
        ans <- rbind(mn, sd)
        dimnames(ans) <- list(c("mean", "sd"), nm)
        class(ans) <- "mcmcMean"
        return(ans)
    } else {
        if ( !v ) names(mn) <- nm
        return(mn)
    }
}

mcmcMeans <- mcmcMean

print.mcmcMean <- function(x, digits = getOption("digits"), ...)
{
    ans <- format(x, digits = 3)
    ans[1, ] <- sapply(ans[1, ], function(x) paste("", x))
    ans[2, ] <- sapply(ans[2, ], function(x) paste("(", x, ")", sep = ""))
    dn <- dimnames(ans)
    dn[[1]] <- rep("", 2)
    if (is.null(dn[[2]]))
        dn[[2]] <- paste('V', 1:dim(ans)[2], sep='.')
    dn[[2]] <- paste(substring("      ", 1, (nchar(ans[2, ]) -
        nchar(dn[[2]]))%/%2), dn[[2]])
    dn[[2]] <- paste(dn[[2]], substring("      ", 1, (nchar(ans[2,
        ]) - nchar(dn[[2]]))%/%2))
    dimnames(ans) <- dn
    print(ans, quote = FALSE)
    return(x)
}

ergMean <- function(x, m = 1)
{
    ### ergodic means
    if (hasArg(m)) m <- max(1, round(m))
    n <- NROW(x)
    if ( m > n )
        stop("Need m <= n")
    if ( is.null(dm <- dim(x)) || dm[2] == 1 )
    {
        ## univariate
        if ( m == 1 )
            ans <- cumsum(x) / 1:n
        else
            if ( m == n )
                ans <- mean(x)
            else
                ans <- cumsum(c(sum(x[1:m]), x[(m+1):n])) / m:n
    } else {
        if ( length(dm) == 2 )
        {
            ## matrix "x" - multivariate
            nm <- colnames(x)
            if ( is.null(nm) )
                nm <- paste(deparse(substitute(x)), 1:NCOL(x), sep='.')
            if ( m == 1 )
                ans <- apply(x, 2, cumsum) / 1:n
            else
                if ( m == n )
                    ans <- matrix(colMeans(x), nrow = 1)
                else
                    ans <- apply(rbind(colSums(x[1:m,]), x[(m+1):n,]), 2, cumsum) / m:n
            colnames(ans) <- nm
        }
        else
            stop("\'x\' must be a vector or a matrix")
    }
    return(ans)
}


###### drop first element of a ts or mts retaining the class
dropFirst <- function(x)
{
    st <- start(x) + c(0, 1)
    freq <- frequency(x)
    newStart <- c(st[1], st[2] %% freq) + c(st[2] %/% freq, 0)
    y <- window(x, start = newStart)
    if (is.null(tsp(x)))
        attr(y, "tsp") <- NULL
    return(y)
}

######
###### Gibbs sampler for "d-inverse-gamma" model
######
dlmGibbsDIG <- function(y, mod, a.y, b.y, a.theta, b.theta, shape.y, rate.y,
                        shape.theta, rate.theta, n.sample = 1,
                        thin = 0, ind, save.states = TRUE,
                        progressBar = interactive())
##################################################################################
##################################################################################
### Gibbs sampler for the 'd-inverse-gamma' model                              ###
### Constant DLMs and univariate observations only                             ###
###                                                                            ###
### y       : data (vector or univariate time series).                         ###
### mod     : a dlm model for the data, with a diagonal 'W' component.         ###
### a.y     : prior mean of the observation precision.                         ###
### b.y     : prior variance of the observation precision.                     ###
### a.theta : vector of prior mean(s) of the system precision(s);              ###
###           recycled if needed.                                              ###
### b.theta : vector of prior variance(s) of the system precision(s);          ###
###           recycled if needed.                                              ###
### shape.y, rate.y : shape and rate parameters of the prior gamma             ###
###           distribution of the observation precision. Can be specified      ###
###           in alternative to 'a' and 'b'.                                   ###
### shape.theta, rate.theta : vectors of shape and rate parameters of the      ###
###           prior distributions of the system precision(s). Can be           ###
###           specified in alternative to 'alpha' and 'beta'.                  ###
### n.sample : number of simulated values in the output.                       ###
### thin    : number of sweeps to discard between each pair of returned        ###
###           draws.                                                           ###
### ind     : vector of indices. If specified, the sampler will only draw      ###
###           the diagonal elements of 'W' having the specified indices.       ###
###           Useful when some of the system variances are zero.               ###
### save.states : if TRUE, the generated states will be returned together      ###
###           with the generated parameter values.                             ###
### progressBar : if TRUE, a text progress bar will be displayed               ###
###           during execution                                                 ###
###                                                                            ###
### Value:                                                                     ###
### A list with components 'dV', 'dW', 'theta' (only if 'save.states' is       ###
### TRUE). 'dV' contains the generated observation variances, 'dW' the         ###
### generated system variances, 'theta' the generated states.                  ###
##################################################################################
##################################################################################
{
    msg1 <- "Either \"a.y\" and \"b.y\" or \"shape.y\" and \"rate.y\" must be specified"
    msg2 <- "Unexpected length of \"shape.y\" and/or \"rate.y\""
    msg3 <- "Unexpected length of \"a.y\" and/or \"b.y\""
    msg4 <- paste("Either \"a.theta\" and \"b.theta\" or \"shape.theta\"",
                  "and \"rate.theta\" must be specified")
    msg5 <- "Unexpected length of \"shape.theta\" and/or \"rate.theta\""
    msg6 <- "Unexpected length of \"a.theta\" and/or \"b.theta\""
    msg7 <- "\"thin\" must be a nonnegative integer"
    msg8 <- "multivariate observations are not allowed"
    msg9 <- "inadmissible value of \"ind\""
    if ( NCOL(y) > 1 )
        stop(msg8)
    r <- ncol(mod$FF)
    if ( hasArg(ind) ) {
        ind <- unique(as.integer(ind))
        s <- 1:r
        if ( !all(ind %in% s) )
            stop(msg9)
        perm <- s[c(ind, s[ !(s %in% ind)])]
        FF(mod) <- mod$FF[, perm, drop = FALSE]
        GG(mod) <- mod$GG[perm, perm, drop = FALSE]
        W(mod) <- mod$W[perm, perm, drop = FALSE]
        p <- length(ind)
    }
    else {
        perm <- ind <- 1 : r
        p <- r
    }
    nobs <- NROW(y)
    if ( is.numeric(thin) && (thin <- as.integer(thin)) >= 0 )
    {
        every <- thin + 1
        mcmc <- n.sample * every
    }
    else
        stop(msg7)
    ## check hyperpriors for precision of 'y'
    if ( !hasArg(a.y) )
        if ( !hasArg(shape.y) ) stop(msg1)
        else
            if ( !hasArg(rate.y) ) stop(msg1)
            else
            {
                ## check length of shape.y and rate.y
                if (!all(c(length(shape.y), length(rate.y)) == 1))
                    warning(msg2)
            }
    else
        if ( !hasArg(b.y) ) stop(msg1)
        else
        {
            if (!all(c(length(a.y), length(b.y)) == 1))
                warning(msg3)
            shape.y <- a.y^2 / b.y
            rate.y <- a.y / b.y
        }
    ## check hyperpriors for precision(s) of 'theta'
    if ( !hasArg(a.theta) )
        if ( !hasArg(shape.theta) ) stop(msg4)
        else
            if ( !hasArg(rate.theta) ) stop(msg4)
            else
            {
                ## check length of shape.theta and rate.theta
                if (!all(c(length(shape.theta), length(rate.theta)) %in% c(1,p)))
                    warning(msg5)
            }
    else
        if ( !hasArg(b.theta) ) stop(msg4)
        else
        {
            if (!all(c(length(a.theta), length(b.theta)) %in% c(1,p)))
                warning(msg6)
            shape.theta <- a.theta^2 / b.theta
            rate.theta <- a.theta / b.theta
        }
    shape.y <- shape.y + 0.5 * nobs
    shape.theta <- shape.theta + 0.5 * nobs
    shape.theta <- rep(shape.theta, length.out = p)
    rate.theta <- rep(rate.theta, length.out = p)
    theta <- matrix(0, nobs + 1, r)
    if ( save.states )
        gibbsTheta <- array(0, dim = c(nobs + 1, r, n.sample))
    gibbsV <- vector("numeric", n.sample)
    gibbsW <- matrix(0, nrow = n.sample, ncol = p)
    it.save <- 0
    if (progressBar) pb <- txtProgressBar(0, mcmc, style = 3)
    for (it in 1:mcmc)
    {
        if (progressBar) setTxtProgressBar(pb, it)
        ## generate states - FFBS
        modFilt <- dlmFilter(y, mod, simplify=TRUE)
        theta[] <- dlmBSample(modFilt)
        ## generate V
        y.center <- y - tcrossprod(theta[-1,,drop=FALSE], mod$FF)
        SSy <- drop(crossprod(y.center))
        mod$V[] <- 1 / rgamma(1, shape = shape.y,
                              rate = rate.y + 0.5 * SSy)
        ## generate W
        theta.center <- theta[-1,,drop=FALSE] -
            tcrossprod(theta[-(nobs + 1),,drop=FALSE], mod$GG)
        SStheta <- drop(sapply( 1 : p, function(i) crossprod(theta.center[,i])))
        SStheta <- colSums((theta[-1,1:p,drop=FALSE] -
                            tcrossprod(theta[-(nobs + 1),,drop=FALSE],mod$GG)[,1:p])^2)
        diag(mod$W)[1:p] <-
            1 / rgamma(p, shape = shape.theta, rate = rate.theta + 0.5 * SStheta)
        ## save
        if ( !(it %% every) )
        {
            it.save <- it.save + 1
            if ( save.states )
                gibbsTheta[,,it.save] <- theta
            gibbsV[it.save] <- diag(mod$V)
            gibbsW[it.save,] <- diag(mod$W)[1:p]
        }
    }
    colnames(gibbsW) <- paste("W", ind, sep = '.')
    if (progressBar) close(pb)
    if ( save.states )
        return(list(dV = gibbsV, dW = gibbsW,
                    theta = gibbsTheta[, order(perm), , drop = FALSE]))
    else
        return(list(dV = gibbsV, dW = gibbsW))
}

