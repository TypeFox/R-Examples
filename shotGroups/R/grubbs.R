#####---------------------------------------------------------------------------
## Grubbs-Patnaik / Grubbs-Pearson / Grubbs-Liu chi^2 distributions
#####---------------------------------------------------------------------------

#####---------------------------------------------------------------------------
## determine Grubbs parameters
#####---------------------------------------------------------------------------

## from eigenvalues of covariance matrix = variance of decorrelated data
## for all dimensionalities
## not vectorized
getGrubbsParam <-
function(sigma, ctr, accuracy=FALSE) {
    sigma <- as.matrix(sigma)
    if(!is.numeric(sigma)) { stop("sigma must be numeric") }
    if(!isTRUE(all.equal(sigma, t(sigma)))) { stop("sigma must be symmetric") }
    if(missing(ctr))       { ctr <- numeric(ncol(sigma)) }

    eig    <- eigen(sigma)               # eigen-system of covariance matrix
    lambda <- eig$values                 # eigenvalue vector
    eVec   <- eig$vectors                # eigenvector matrix

    if(!all(lambda >= -sqrt(.Machine$double.eps) * abs(lambda[1]))) {
        stop("sigma is numerically not positive definite")
    }

    ## take systematic location bias into account?
    ## if so -> rotate ctr with eigenvectors (decorrelated data)
    ## non-centrality parameters delta = ai^2
    deltaI <- if(accuracy) {
        (t(eVec) %*% ctr)^2 / lambda
    } else {
        0
    }

    ## distribution parameters
    m <- sum(lambda) + sum(lambda*deltaI)          # mean (c1)
    v <- 2*sum(lambda^2) + 4*sum(lambda^2*deltaI)  # variance (2*c2)
    n <- 2*m^2 / v                                 # df chi-square Patnaik

    ## for Pearson approximation
    mu3    <- 8*sum(lambda^3 + 3*lambda^3*deltaI)  # 3rd moment = skewness (8*c3)
    beta1  <- mu3^2 / v^3
    nPrime <- 8/beta1                              # df chi-square Pearson (c2^3/c3^2)

    ## for Liu, Tang & Zhang approximation
    c1 <- m
    c2 <- v/2
    c3 <- sum(lambda^3) + 3*sum(lambda^3*deltaI)
    c4 <- sum(lambda^4) + 4*sum(lambda^4*deltaI)
    s1 <- c3/sqrt(c2^3)
    s2 <- c4/c2^2

    ## if(s1^2 > s2)
    a1     <- suppressWarnings(1 / (s1 - sqrt(s1^2 - s2))) # NaN if s1^2 < s2
    delta1 <- s1*a1^3 - a1^2
    l1     <- a1^2 - 2*delta1

    ## if(s1^2 < s2)
    a2     <- 1/s1
    delta2 <- 0
    l2     <- c2^3/c3^2

    a     <- ifelse(s1^2 > s2, a1,     a2)
    delta <- ifelse(s1^2 > s2, delta1, delta2)
    l     <- ifelse(s1^2 > s2, l1    , l2)
    muX   <- l + delta
    varX  <- 2*a^2

    return(list(m=m, v=v, n=n, nPrime=nPrime, muX=muX, varX=varX, l=l, delta=delta))
}

## determine Grubbs parameters (accuracy=FALSE) from Hoyt parameters - just 2D
## vectorized
getGPfromHP <-
function(qpar, omega) {
    nnaQ <- which(!is.na(qpar))
    nnaO <- which(!is.na(omega))
    stopifnot(all(qpar[nnaQ] > 0), all(qpar[nnaQ] < 1), all(omega[nnaO] > 0))

    ## eigenvalues from q and omega
    ev2 <- omega / ((1/qpar^2) + 1)
    ev1 <- omega - ev2

    ## distribution parameters
    m <- ev1+ev2                         # mean
    v <- 2*(ev1^2 + ev2^2)               # variance
    n <- 2*m^2 / v                       # df chi-square Patnaik

    ## for Pearson approximation
    mu3    <- 8*(ev1^3 + ev2^3)          # 3rd moment
    beta1  <- mu3^2 / v^3
    nPrime <- 8/beta1                    # df chi-square Pearson

    ## for Liu, Tang & Zhang approximation
    c1 <- m
    c2 <- v/2
    c3 <- ev1^3 + ev2^3
    c4 <- ev1^4 + ev2^4
    s1 <- c3/sqrt(c2^3)
    s2 <- c4/c2^2

    ## if(s1^2 > s2)
    a1     <- suppressWarnings(1 / (s1 - sqrt(s1^2 - s2))) # NaN if s1^2 < s2
    delta1 <- s1*a1^3 - a1^2
    l1     <- a1^2 - 2*delta1

    ## if(s1^2 < s2)
    a2     <- 1/s1
    delta2 <- 0
    l2     <- c2^3/c3^2

    a     <- ifelse(s1^2 > s2, a1,     a2)
    delta <- ifelse(s1^2 > s2, delta1, delta2)
    l     <- ifelse(s1^2 > s2, l1    , l2)
    muX   <- l + delta
    varX  <- 2*a^2

    return(list(m=m, v=v, n=n, nPrime=nPrime, muX=muX, varX=varX, l=l, delta=delta))
}

## determine Grubbs parameters (accuracy=TRUE) from Rice parameters - just 2D
## vectorized
getGPfromRP <-
function(nu, sigma) {
    nnaN <- which(!is.na(nu))
    nnaS <- which(!is.na(sigma))
    stopifnot(all(nu[nnaN] >= 0), all(sigma[nnaS] > 0))

    ## eigenvalues from sigma
    ev1 <- sigma^2
    ev2 <- sigma^2

    ## take systematic location bias into account
    ## non-centrality parameters delta = ai^2
    ## distribution is circular -> offset in only 1 direction
    deltaI1 <- nu^2/ev1
    deltaI2 <- rep(0, max(c(length(nu), length(sigma))))

    ## distribution parameters
    m <- ev1+ev2 + ev1*deltaI1 + ev2*deltaI2       # mean
    v <- 2*(ev1^2 + ev2^2) + 4*(ev1^2*deltaI1 + ev2^2*deltaI2) # variance
    n <- 2*m^2 / v                                 # df chi-square Patnaik

    ## for Pearson approximation
    ## 3rd moment = skewness
    mu3    <- 8*(ev1^3+ev2^3 + 3*ev1^3*deltaI1 + 3*ev2^3*deltaI2)
    beta1  <- mu3^2 / v^3
    nPrime <- 8/beta1                              # df chi-square Pearson

    ## for Liu, Tang & Zhang approximation
    c1 <- m
    c2 <- v/2
    c3 <- ev1^3+ev2^3 + 3*ev1^3*deltaI1 + 3*ev2^3*deltaI2
    c4 <- ev1^4+ev2^4 + 4*ev1^4*deltaI1 + 4*ev2^4*deltaI2
    s1 <- c3/sqrt(c2^3)
    s2 <- c4/c2^2

    ## if(s1^2 > s2)
    a1     <- suppressWarnings(1 / (s1 - sqrt(s1^2 - s2))) # NaN if s1^2 < s2
    delta1 <- s1*a1^3 - a1^2
    l1     <- a1^2 - 2*delta1

    ## if(s1^2 < s2)
    a2     <- 1/s1
    delta2 <- 0
    l2     <- c2^3/c3^2

    a     <- ifelse(s1^2 > s2, a1,     a2)
    delta <- ifelse(s1^2 > s2, delta1, delta2)
    l     <- ifelse(s1^2 > s2, l1    , l2)
    muX   <- l + delta
    varX  <- 2*a^2

    return(list(m=m, v=v, n=n, nPrime=nPrime, muX=muX, varX=varX, l=l, delta=delta))
}

#####---------------------------------------------------------------------------
## density, cdf, quantile function, and random numbers
#####---------------------------------------------------------------------------

## cdf - vectorized
## change of variables
## Y = u(X)
## v(u(X)) = X -> v = u^-1
## f_Y(y) =  f_X(v(y)) * v'(y), u increasing
## f_Y(y) = -f_X(v(y)) * v'(y), u decreasing
dChisqGrubbs <-
function(x, m, v, n, nPrime, muX, varX, l, delta,
    type=c("Patnaik", "Pearson", "Liu")) {
    type <- match.arg(type)

    dens <- if(type == "Patnaik") {
        uInv  <- x^2*2*m/v
        uInvD <- x*4*m/v                # derivative of u^-1
        dchisq(uInv, df=n) * uInvD
    } else if(type == "Pearson") {
        uInv  <- (x^2-m)*sqrt(2*nPrime/v) + nPrime
        uInvD <- x*2*sqrt(2*nPrime/v)
        dchisq(uInv, df=nPrime) * uInvD
    } else if(type == "Liu") {
        uInv  <- (x^2-m)*sqrt(varX/v) + muX
        uInvD <- x*2*sqrt(varX/v)
        dchisq(uInv, df=l, ncp=delta) * uInvD
    }

    return(dens)
}

## cdf - vectorized
pChisqGrubbs <-
function(q, m, v, n, nPrime, muX, varX, l, delta, lower.tail=TRUE,
    type=c("Patnaik", "Pearson", "Liu")) {
    type <- match.arg(type)

    pp <- if(type == "Patnaik") {
        qq <- q^2*2*m/v
        pchisq(qq, df=n, lower.tail=lower.tail)
    } else if(type == "Pearson") {
        qq <- (q^2-m)*sqrt(2*nPrime/v) + nPrime
        pchisq(qq, df=nPrime, lower.tail=lower.tail)
    } else if(type == "Liu") {
        qq <- (q^2-m)*sqrt(varX/v) + muX
        pchisq(qq, df=l, ncp=delta, lower.tail=lower.tail)
    }

    return(pp)
}

## quantile function - vectorized
qChisqGrubbs <-
function(p, m, v, n, nPrime, muX, varX, l, delta, lower.tail=TRUE,
    type=c("Patnaik", "Pearson", "Liu")) {
    type <- match.arg(type)

    qq <- if(type == "Patnaik") {
        q1 <- qchisq(p, df=n,            lower.tail=lower.tail)
        sqrt((v/(2*m)) * q1)
    } else if(type == "Pearson") {
        q1 <- qchisq(p, df=nPrime,       lower.tail=lower.tail)
        sqrt((q1-nPrime)*sqrt(v/(2*nPrime)) + m)
    } else if(type == "Liu") {
        q1 <- qchisq(p, df=l, ncp=delta, lower.tail=lower.tail)
        sqrt((q1-muX) * sqrt(v/varX) + m)
    }

    return(qq)
}

## random numbers from Grubbs chi^2 distributions
rChisqGrubbs <-
function(nr, m, v, n, nPrime, muX, varX, l, delta,
         type=c("Patnaik", "Pearson", "Liu")) {
    u <- runif(nr)                       # uniform random numbers
    qChisqGrubbs(u, m=m, v=v, n=n, nPrime=nPrime, muX=muX,
                 varX=varX, l=l, delta=delta, type=type)
}
