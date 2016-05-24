Wilks.test <- function(x, ...) UseMethod("Wilks.test")

Wilks.test.formula <- function(formula, data, ..., subset, na.action)
{
    m <- match.call(expand.dots = FALSE)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    if (.check_vars_numeric(m))
        stop("Wilks test applies only to numerical variables")

    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)

    xint <- match("(Intercept)", colnames(x), nomatch=0)
    if(xint > 0) x <- x[, -xint, drop=FALSE]

    res <- Wilks.test.default(x, grouping, ...)
    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1]] <- as.name("Wilks.test")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- .getXlevels(Terms, m)
    res$na.action <- attr(m, "na.action")
    res
}

Wilks.test.data.frame <- function(x, ...)
{
    res <- Wilks.test(structure(data.matrix(x), class="matrix"), ...)
    cl <- match.call()
    cl[[1]] <- as.name("Wilks.test")
    res$call <- cl
    res
}

Wilks.test.matrix <- function(x, grouping, ..., subset, na.action)
{
    if(!missing(subset)) {
        x <- x[subset, , drop = FALSE]
        grouping <- grouping[subset]
    }
    if(!missing(na.action)) {
        dfr <- na.action(structure(list(g = grouping, x = x),
                                   class = "data.frame"))
        grouping <- dfr$g
        x <- dfr$x
    }

    res <- Wilks.test.default(x, grouping, ...)
    cl <- match.call()
    cl[[1]] <- as.name("Wilks.test")
    res$call <- cl
    res
}

##
## Default S3 method for Wilks.test
##
##  x   - matrix or data frame
##  grp - grouping variable - factor specifying the class for each observation
##  approximation   -
##  method  - "c" for standard estimators of the mean and variance, "mcd" for robust estimates based on MCD.
##  xq  - used only when method="mcd". Multiplication factor for the approximate Chi2 distribution. If NULL
##          simulation will be performed to determine xq and xd (default is NULL)
##  xd  - used only when method="mcd". Degrees of freedom for the approximate Chi2 distribution. If NULL
##          simulation will be performed to determine xq and xd (default is NULL)
##  nrep - number of trials for the simulation for estimating xq and xd (default is 3000)
##  trace   -   whether to provide trace output (default is FALSE)
##
Wilks.test.default <- function(x,
                            grouping,
                            method=c("c", "mcd", "rank"),
                            approximation=c("Bartlett", "Rao", "empirical"), 
                            xd=NULL, xq=NULL, xfn = NULL, xwl=NULL, 
                            nrep=3000,
                            trace=FALSE, ...){

    alpha <- 0.5
##    approximation <- "Bartlett"
##    approximation <- "empirical"

    cl <- match.call()
    method <- match.arg(method)
    approximation <- match.arg(approximation)
    dname <- deparse(substitute(x))

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
        stop("All observations have missing values!")

    grouping <- grouping[ok]

    n <- nrow(x)
    p <- ncol(x)

    if(method == "rank")
        x <- apply(x, 2, rank)

    ww <- .wilks(x, grouping, method, alpha=alpha)
    nk <- ww$nk
    ng <- length(nk)
    wilks <- ww$wilks

    METHOD = NULL
    PARAMETER = NULL
    if(approximation == "Bartlett") {

        ## Bartlett's ChiSq approximation
        Y <- log(wilks)
        if(method == "mcd"){
            if(missing(xd) || is.null(xd) || missing(xq) || is.null(xq)){
                xqd <- simulateChi2(nrep=nrep, nk, p, trace=trace, alpha=alpha)
                xd <- xqd$xd
                xq <- xqd$xq
                xfn <- xqd$xfn
                xwl <- xqd$xwl
            }
            STATISTIC <- wilks
            chi2.value <- Y/xd
            df <- xq
            names(STATISTIC) <- "Wilks' Lambda"     # "Chi2-approximation (simulated)"
            METHOD <- "Robust One-way MANOVA (Bartlett Chi2)"
        }else{
            STATISTIC <- wilks
            chi2.value <-  -(n - 1 - (p+ng)/2)*Y
            df <- p*(ng-1)
            names(STATISTIC) <- "Wilks' Lambda"     # "Chi2-approximation"
            METHOD <- "One-way MANOVA (Bartlett Chi2)"
        }
        PVAL <- 1-pchisq(chi2.value, df)
        PARAMETER <- c(chi2.value, df)
        names(PARAMETER) <- c("Chi2-Value", "DF")
    }else if(approximation == "Rao"){
## Rao's F approximation
        t1 <- p*p*(ng-1)*(ng-1) - 4
        t2 <- p*p + (ng-1)*(ng-1) - 5
        t <- if(t1==0 | t2 == 0) 1 else sqrt(t1/t2)
##      t <- sqrt((p*p*(ng-1)*(ng-1) - 4)/(p*p + (ng-1)*(ng-1) - 5))

        df1 <- p*(ng-1)
        df2 <- t*((n-1) - (p + ng)/2) - (p*(ng-1) - 2)/2        #   instead of n - the total number of observations - use
                                                                #   the sum of weights ?
        Y <- wilks^(1/t)
        STATISTIC <- df2*(1-Y)/(df1*Y)
        PVAL <- 1-pf(STATISTIC, df1, df2)
        PARAMETER <- c(df1,df2)
        names(STATISTIC) <- "F-approximation"
        names(PARAMETER) <- c("DF1", "DF2")
    }else if(approximation == "empirical") {
        Y <- log(wilks)
        if(method == "mcd"){
            if(missing(xd) || is.null(xd) || missing(xq) || is.null(xq)){
                xqd <- simulateChi2(nrep=nrep, nk, p, trace=trace, alpha=alpha)
                xd <- xqd$xd
                xq <- xqd$xq
                xfn <- xqd$xfn
                xwl <- xqd$xwl
            }
            STATISTIC <- wilks
            METHOD <- "Robust One-way MANOVA (empirical distribution)"
            PVAL <- xfn(wilks)

        }else{
            STATISTIC <- wilks
            chi2.value <-  -(n - 1 - (p+ng)/2)*Y
            df <- p*(ng-1)
            names(STATISTIC) <- "Wilks' Lambda"     # "Chi2-approximation"
            METHOD <- "One-way MANOVA (Bartlett Chi2)"
            PVAL <- 1-pchisq(chi2.value, df)
            PARAMETER <- c(chi2.value, df)
            names(PARAMETER) <- c("Chi2-Value", "DF")
        }
    }

    ans <- list(statistic=STATISTIC,
                parameter=PARAMETER,
                p.value=PVAL,
                estimate=ww$group.means,
                method=METHOD,
                data.name=dname,
                W=ww$W,
                T=ww$T,
                wilks=wilks,
                xd=xd,
                xq=xq,
                xfn=xfn,
                xwl=xwl)

    cl[[1]] <- as.name("Wilks.test")
    ans$call <- cl
    class(ans) <- "htest"
    ans
}

covMWcd <-function(x, grouping, alpha=0.75){

    ## compute group means, pool the observations and compute the common
    ##  covariance matrix
    covMWcd.B <- function(){
        group.means <- matrix(0, ng, p)
        for(i in 1:ng){
            mcd <- CovMcd(x[which(grouping == lev[i]),], alpha=alpha)
            group.means[i,] <- getCenter(mcd)
        }

        mcd <- CovMcd(x - group.means[g,], alpha=alpha)

        ans <- list(center=group.means, wcov=getCov(mcd), mah=getDistance(mcd), wt=mcd@wt)

        class(ans) <- "mwcd"
        attr(ans, "call") <- sys.call()
        return(ans)
    }

    x <-as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)

    g <- as.factor(grouping)
    lev <- levels(g)
    ng <- length(lev)

    covMWcd.B()
 }

##
##  Find by simulation the parameters 'd' and 'q' of the approximate Chi2 distribution
##  for the robust Wilks Lambda test based on MCD
##  for a given p, g and n=\sum ni
##
##  nrep    - number of trials (3000)
##  nk      - array of integers giving the group sizes, e.g. nk=c(20, 20)
##  p       - dimension
##
simulateChi2 <- function(nrep=3000, nk, p, trace=FALSE, alpha=0.75){
    if(missing(nk))
        stop("The number and size of groups must be provided!")
    if(missing (p))
        stop("The dimensipon 'p' must be provided")

    if(trace)
        cat("\nFind approximate distribution...\n")

    ptm <- proc.time()

    n <- sum(nk)
    ng <- length(nk)
    grp <- as.factor(rep(1:ng, nk))
    wl <- rwl <- vector(mode = "numeric", length = nrep)

    for(i in 1:nrep){
#        x <- as.matrix(mvrnorm(n=n, mu=rep(0,p), Sigma=diag(rep(1,p))))
        x <- matrix(rnorm(n*p), ncol=p)

        tt <- .wilks(x, grp, method="c", alpha=alpha)
        rtt <- .wilks(x, grp, method="mcd", alpha=alpha)

        wl[i] <- tt$wilks
        rwl[i] <- rtt$wilks
    }

    Y <- log(wl)
    RY <- log(rwl)

    xmean <- mean(Y)
    xvar <- var(Y)
    xq <- 2*xmean*xmean/xvar
    xd <- xmean/xq

    rxmean <- mean(RY)
    rxvar <- var(RY)
    rxq <- 2*rxmean*rxmean/rxvar
    rxd <- rxmean/rxq
    rxfn <- ecdf(rwl)

    if(trace){
        cat(sprintf('\n  nrep=%5.0f p=%4.0f ng=%2.0f n=%5.0f', nrep,p,ng,n))
        cat(sprintf('\n  alpha=%5.2f', alpha))
        cat(sprintf('\n         mean   var    d     q'))
        cat(sprintf('\n                      %5.3f %5.3f', -1/(n-1-(p*ng)/2), p*(ng-1)))
        cat(sprintf('\n  WILKS: %5.3f %5.3f %5.3f %5.3f', xmean, xvar, xd, xq))
        cat(sprintf('\nR WILKS: %5.3f %5.3f %5.3f %5.3f', rxmean, rxvar, rxd, rxq))
        cat("\n")
        cat(" Elapsed time: ", (proc.time() - ptm)[1],"\n")
    }
    invisible(list(xd=rxd, xq=rxq, xfn=rxfn, xwl=rwl))
}

.wilks <- function(x, grouping, method, alpha=0.75){
    n <- nrow(x)
    p <- ncol(x)
    if(!is.factor(g <- grouping))
        g <- as.factor(grouping)

    ## get rid of empty groups
    grouping <- g <- as.factor(as.character(g))

    lev <- levels(g)
    counts <- as.vector(table(g))
    if(any(counts == 0)) {
        stop(paste("group(s)", paste(lev[counts == 0], collapse=" "),"are empty"))
    }
    ng <- length(counts)

    if(method == "c" | method == "rank"){
        wts <- rep(1, n)
    }else if(method == "mcd"){
        mcd <- covMWcd(x, grouping=grouping, alpha=alpha)
        wts <- mcd$wt
    }else
        stop("Undefined method: ", method)

    group.means <- matrix(0,ng,p)
    for(i in 1:ng){
        group.means[i,] <- cov.wt(x[which(grouping == lev[i]),], wt=wts[which(grouping == lev[i])])$center
    }

    wcross <- cov.wt((x - group.means[g, ]), wt=wts)
    wcross <- (sum(wts)-1) * wcross$cov
    tcross <- cov.wt(x, wt=wts)
    tcross <- (sum(wts)-1) * tcross$cov

    wilks <- det(wcross)/det(tcross)
    names(wilks) <- "Wilks' Lambda"
    names(group.means) <- names(x)

    dimnames(group.means) <- list(lev, dimnames(x)[[2]])

    list(nk=counts,
         wilks=wilks,
         W=wcross,
         T=tcross,
         group.means=group.means)
}

.wilksx <- function(x, grouping, method, alpha=0.75){
    n <- nrow(x)
    p <- ncol(x)
    if(!is.factor(g <- grouping))
        g <- as.factor(grouping)

    ## get rid of empty groups
    g <- as.factor(as.character(g))
    lev <- levels(g)
    counts <- as.vector(table(g))
    if(any(counts == 0)) {
        stop(paste("group(s)", paste(lev[counts == 0], collapse=" "),"are empty"))
    }
    ng <- length(counts)

    if(method == "c" | method == "rank"){
        wts <- rep(1, n)
    }else if(method == "mcd"){
        mcd <- covMWcd(x, grouping=grouping, alpha=alpha)
        wts <- mcd$wt
    }else
        stop("Undefined method: ", method)

    group.means <- matrix(0,ng,p)
    for(i in 1:ng){
        group.means[i,] <- cov.wt(x[which(grouping == lev[i]),], wt=wts[which(grouping == lev[i])])$center
    }

    wcross <-cov.wt((x - group.means[g, ]), wt=wts)
    wcross <- (sum(wts)-1) * wcross$cov
    tcross <- cov.wt(x, wt=wts)
    tcross <- (sum(wts)-1) * tcross$cov

    wilks <- det(wcross)/det(tcross)
    names(wilks) <- "Wilks' Lambda"
    names(group.means) <- names(x)

    dimnames(group.means) <- list(lev, dimnames(x)[[2]])

    list(nk=counts,
         wilks=wilks,
         W=wcross,
         T=tcross,
         group.means=group.means)
}
