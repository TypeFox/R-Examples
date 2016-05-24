T2.test <- function(x, ...) UseMethod("T2.test")

T2.test.default <- function(x, y = NULL, mu = 0, conf.level = 0.95, method=c("c", "mcd"), ...)
{
    xcall <- match.call()
    method <- match.arg(method)

    if(!is.null(y)){
        dname <- paste(deparse(substitute(x)),"and",
                        deparse(substitute(y)))
    } else
        dname <- deparse(substitute(x))

    ##  Validate that x is a numerical matrix or a dataframe,
    ##  convert to a matrix and drop all rows with missing values
    x <- .tomatrix(x)
    dx <- dim(x)
    n <- n1 <- dx[1]
    p <- dx[2]

    if(!is.null(y))
    {

        ##  Validate that y is a numerical matrix or a dataframe,
        ##  convert to a matrix and drop all rows with missing values
        y <- .tomatrix(y)
        dy <- dim(y)
        n2 <- dy[1]
        py <- dy[2]
        if(p != py)
            stop("'x' and 'y' must have the same dimension!")
    }


    if(!is.numeric(mu) || ((lmu <- length(mu)) > 1 & lmu != p))
        stop("'mu' must be a numeric vector of length ", p)
    if(lmu == 1)
        mu <- rep(mu, p)

    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")
    alpha <- 1 - conf.level

    ## One sample test -  y == NULL
    if(is.null(y)) {
        alternative <- paste("true mean vector is not equal to",
                    paste("(", paste(round(mu, digits=3), collapse = ", "), ")'", sep = ""), "\n")
        null.value <- mu
        METHOD = "One-sample Hotelling test"

        if(method == "c"){
            xbar <- colMeans(x)
            xdiff <- xbar - mu
            V <- var(x)

            d <- p * (n-1)/(n-p)
            q <- n-p
            T2 <- STATISTIC <- n * crossprod(xdiff, solve(V, xdiff))[1, ]
            F  <- STATISTIC <- STATISTIC/d

            PVALUE <- 1 - pf(STATISTIC, p, n - p)
            PARAMETER <- c(p, n - p)
            ESTIMATE = t(as.matrix(xbar))
            rownames(ESTIMATE) <- "mean x-vector"
        }else if(method == "mcd"){
            mcd <- CovMcd(x, trace=FALSE, alpha=0.75)
            xbar <- getCenter(mcd)
            xdiff <- xbar - mu
            V <- getCov(mcd)

            xdq <- .getApprox(p, n)
            d <- xdq$d
            q <- xdq$q

            T2 <- STATISTIC <- n * crossprod(xdiff, solve(V, xdiff))[1, ]
            F  <- STATISTIC <- STATISTIC/d

            PVALUE <- 1 - pf(STATISTIC, p, q)
            PARAMETER <- c(p, q)
            ESTIMATE = t(as.matrix(xbar))
            rownames(ESTIMATE) <- "MCD x-vector"
            METHOD <- paste(METHOD, " (Reweighted MCD Location)")
        } else
            stop(paste("Invalid method=",method))

        cutoff.alpha <- qf(1-alpha, p, q)
        ## simultaneous confidence intervals for the components of mu
        conf.int <- matrix(1:(2*p), nrow=p)
        for(i in 1:p) {
            conf.int[i,1] <- xbar[i] - sqrt(1/n*d*qf(1-alpha, p, q) * V[i,i])
            conf.int[i,2] <- xbar[i] + sqrt(1/n*d*qf(1-alpha, p, q) * V[i,i])
        }
        dimnames(conf.int) <- list(dimnames(x)[[2]], c("Lower bound","Upper bound"))
        attr(conf.int,"conf.level") <- conf.level

        ## switch off the confidence intervals, since 'htest'
        ## does not know how to print them
        conf.int <- NULL
    } else {
        if(method != "c")
            stop("Robust two-sample test not yet implemeted!")

        xbar <- colMeans(x)
        ybar <- colMeans(y)
        xdiff <- xbar - ybar                                    # the difference between the two means
        Vx <- var(x)
        Vy <- var(y)
        V <- ((n1 - 1) * Vx + (n2 - 1) * Vy) / (n1+ n2 - 2)     # the pooled covariance matrix

        df1 <- p
        df2 <- n1 + n2 - p - 1
        T2 <- STATISTIC <- crossprod(xdiff, solve(V, xdiff))[1,] * n1 * n2 / (n1+n2)
        F  <- STATISTIC <- STATISTIC * (n1 + n2 - p - 1) / (n1 + n2 - 2) / p
        PVALUE <- 1 - pf(STATISTIC, df1, df2)
        PARAMETER = c(df1, df2)

        null.value <- NULL
        METHOD = "Two-sample Hotelling test"
        ESTIMATE = rbind(xbar, ybar)
        rownames(ESTIMATE) <- c("mean x-vector", "mean y-vector")
        alternative <- paste("true difference in mean vectors is not equal to (", paste(rep(0,p), collapse=","),")", sep="")

        conf.int <- NULL
    }

    names(PARAMETER) <- c("df1", "df2")
##    names(STATISTIC) <- "T^2"
    STATISTIC <- c(T2, F)
    names(STATISTIC) <- c("T2", "F")


    rval <- list(statistic = STATISTIC,
            parameter = PARAMETER,
            p.value = PVALUE,
            conf.int=conf.int,
            estimate=ESTIMATE,
            null.value = NULL,
            alternative = alternative,
            method=METHOD,
            data.name=dname)

    class(rval) <- "htest"
    return(rval)
}

T2.test.formula <- function(formula, data, subset, na.action, ...)
{
    if(missing(formula)
       || (length(formula) != 3)
       || (length(attr(terms(formula[-2]), "term.labels")) != 1))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)

    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)

    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2)
        stop("grouping factor must have exactly 2 levels")

    xind <- which(g==levels(g)[1])
    yind <- which(g==levels(g)[2])

    y <- T2.test(x=mf[[response]][xind,], y=mf[[response]][yind,], ...)
    y$data.name <- DNAME
    y
}

##
## Validate that 'x' is a numeric data frame or matrix.
## Convert 'x' to a matrix
## Optionally drop all rows with missing values
##
.tomatrix <- function(x, drop.missing=TRUE){
    x.name <- deparse(substitute(x))
    msg <- paste(x.name, " is not a numeric dataframe or matrix.")
    if(is.vector(x) || (is.matrix(x) && !is.data.frame(x))) {
        if(!is.numeric(x))
        stop(msg)
    }
    if(!is.vector(x) && !is.matrix(x) || is.data.frame(x)) {
        if((!is.data.frame(x) && !is.numeric(x)) ||
        (!all(sapply(x, data.class) == "numeric")))
        stop(msg)
    }
    x <- as.matrix(x)
    if(drop.missing){
        ## drop all rows with missing values
        na.x <- !is.finite(x %*% rep(1, ncol(x)))
        xok <- !na.x
        x <- x[xok,  , drop = FALSE]
        dx <- dim(x)
        if(!length(dx))
            stop("All observations have missing values!")
    }
    x
}


.kappa.s <- function(alfa, p)
{
    ## Computes the asymptotic variance for the reweighted
    ##  MCD location estimator
    ##  - p is the dimension
    ##  - alfa, between 1/2 and 1 is the percentage of the observations
    ##      that are used to determine the raw MCD.

    alfa <- 1 - alfa
    qalfa <- qchisq(1 - alfa, p)
    calfainvers <- pgamma(qalfa/2, p/2 + 1)
    c2 <- -1/2 * pgamma(qalfa/2, p/2 + 1)
    effrr <- (4 * c2^2)/calfainvers
    delta <- 0.025
    qdelta <- qchisq(1 - delta, p)
    c1invers <- pgamma(qdelta/2, p/2 + 1)
    c2delta <- -1/2 * pgamma(qdelta/2, p/2 + 1)
    c1tilde <- pgamma(qdelta/2, p/2)
    asvar <- (1 + (2 * c2delta)/c1tilde)^2/effrr
    asvar <- asvar - (calfainvers * (c1tilde + 2 * c2delta))/(c1tilde^2 * c2)
    asvar <- asvar + c1invers/c1tilde^2
    return(asvar)
}

##
##  F approximation for the robust Hotelling test - Willems et al
##  Calculate the factor d and the degrees of freedom q for a given n and p
##
.getApprox <- function(p, n)
{
    ## determining the estimates for E[X] and Var[X]
    coeffpqp2varx <- matrix(c(7,4.829586,2.868279,10,4.261745,3.023672),ncol=2)
    coeffpqp2ex <- matrix(c(7,1.050505,1.360808,10,0.7280273,1.4685161),ncol=2)
    limvarx <- .kappa.s(0.75, p)^2*2*p
    limex <- .kappa.s(0.75, p)*p

    vb1 <- exp(coeffpqp2varx[2,1])/p^(coeffpqp2varx[3,1])
    vb2 <- exp(coeffpqp2varx[2,2])/p^(coeffpqp2varx[3,2])
    vb <- c(log(vb1),log(vb2))
    va12 <- log(7*p^2)
    va22 <- log(10*p^2)
    vA <- matrix(c(1,va12,1,va22), ncol=2, byrow=TRUE)
    vy <-   if(p == 2)          c(8.518773,-1.851881)
            else if(p == 1)     c(6.202284,-1.731468)
            else                solve(vA,vb)
    varx <- limvarx+exp(vy[1])/(n^(-vy[2]))

    eb1 <- exp(coeffpqp2ex[2,1])/p^(coeffpqp2ex[3,1])
    eb2 <- exp(coeffpqp2ex[2,2])/p^(coeffpqp2ex[3,2])
    eb <- c(log(eb1),log(eb2))
    ea12 <- log(7*p^2)
    ea22 <- log(10*p^2)
    eA <- matrix(c(1,ea12,1,ea22), ncol=2, byrow=TRUE)
    ey <- if(p == 2)        c(3.383723,-1.081893)
            else if(p == 1)   c(2.741737,-1.403526)
            else              solve(eA, eb)
    ex <- limex + exp(ey[1])/(n^(-ey[2]))

    ## matching the moments
    q <- (varx/ex^2*p/2-1)^(-1)*(p+2) + 4

    ## When n gets large, the expression for q goes to infinity,
    ## but is very sensitive;  no harm is done by setting q equal
    ## to n when the expression yields negative or extremely large values
    if(q > n || q < 0)
        q <- n

    d <- ex*(q-2)/q
    list(d=d, q=q)
}
