# $Id: wilcox.exact.R,v 1.17 2003/01/16 10:44:25 hothorn Exp $

wilcox.exact <- function(x, ...) UseMethod("wilcox.exact")

wilcox.exact.default <-
function(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
         mu = 0, paired = FALSE, exact = NULL, 
         conf.int = FALSE, conf.level = 0.95, ...) 
{
    alternative <- match.arg(alternative)
    if(!missing(mu) && ((length(mu) > 1) || !is.finite(mu))) 
        stop("mu must be a single number")
    if(conf.int) {
        if(!((length(conf.level) == 1)
             && is.finite(conf.level)
             && (conf.level > 0)
             && (conf.level < 1)))
            stop("conf.level must be a single number between 0 and 1")
    }
    MIDP <- NULL

    if(!is.numeric(x)) stop("`x' must be numeric")

    if(!is.null(y)) {
        if(!is.numeric(y)) stop("`y' must be numeric")
        DNAME <- paste(deparse(substitute(x)), "and",
                       deparse(substitute(y)))
        if(paired) {
            if(length(x) != length(y)) 
                stop("x and y must have the same length")
            OK <- complete.cases(x, y)
            x <- x[OK] - y[OK]
            y <- NULL
        }
        else {
            x <- x[is.finite(x)]
            y <- y[is.finite(y)]
        }
    } else {
        DNAME <- deparse(substitute(x))
        if(paired) 
            stop("y missing for paired test")
        x <- x[is.finite(x)]
    }

    if(length(x) < 1) 
        stop("not enough (finite) x observations")
    CORRECTION <- 0
    if(is.null(y)) {
        METHOD <- "Exact Wilcoxon signed rank test"
        x <- x - mu
        ZEROES <- any(x == 0)
        if(ZEROES) 
            x <- x[x != 0]
        n <- length(x)
        if(is.null(exact)) 
            exact <- (n < 50)
        r <- rank(abs(x))
        STATISTIC <- sum(r[x > 0])
        names(STATISTIC) <- "V"
        TIES <- (length(r) != length(unique(r)))

        if (exact) {
            PVAL <-
               switch(alternative,
                      "two.sided" = pperm(STATISTIC, r, n,
                                          alternative="two.sided", pprob=TRUE),
                      "greater" = pperm(STATISTIC, r, n,
                                        alternative="greater", pprob=TRUE),
                      "less" = pperm(STATISTIC, r, n,
                                     alternative="less", pprob=TRUE))
            MIDP <- PVAL$PPROB
            PVAL <- PVAL$PVALUE
            if(conf.int && !is.na(x)) {
                ## Exact confidence intervale for the median in the
                ## one-sample case.  When used with paired values this
                ## gives a confidence interval for mean(x) - mean(y).
                x <- x + mu             # we want a conf.int for the median
                alpha <- 1 - conf.level
                diffs <- outer(x, x, "+")
                diffs <- sort(diffs[!lower.tri(diffs)]) / 2
                if (TIES) {
                  fs <- function(d) {
                    xx <- x - d; sum(rank(abs(xx))[xx > 0])
                  }
                  w <- sapply(diffs, fs)
                } else {
                  w <- sum(rank(abs(x))):1
                }
                cint <-
                    switch(alternative,
                           "two.sided" = {
                               qu <- qperm(alpha/2, r, n) 
                               ql <- qperm(1-alpha/2, r, n)
                               if (qu <= min(w)) lci <- max(diffs)
                               else lci <- min(diffs[w <= qu])
                               if (ql >= max(w)) uci <- min(diffs)
                               else uci <- max(diffs[w > ql])
                               c(uci, lci)
                           },
                           "greater"= {
                               ql <- qperm(1-alpha, r, n)
                               if (ql >= max(w)) uci <- min(diffs)
                               else uci <- max(diffs[w > ql])
                               c(uci, Inf)
                           },
                           "less"= {
                               qu <- qperm(alpha, r, n)
                               if (qu <= min(w)) lci <- max(diffs)
                               else lci <- min(diffs[w <= qu])
                               c(-Inf, lci)
                })
                attr(cint, "conf.level") <- conf.level    
                wmean <- sum(r)/2
                ESTIMATE <- mean(c(min(diffs[w <= ceiling(wmean)]),
                                   max(diffs[w > wmean])))
                names(ESTIMATE) <- "(pseudo)median"
                }
            } else {
                METHOD <- "Asymptotic Wilcoxon signed rank test"
                wmean <- sum(r)/2
                wvar <- sum(r^2)/4
                PVAL <- pnorm((STATISTIC - wmean) / sqrt(wvar))
                if(alternative == "greater")
                    PVAL <- 1 - PVAL
                if(alternative == "two.sided") 
                    PVAL <- 2 * min(PVAL, 1-PVAL) 

                if(conf.int && !is.na(x)) {
                    ## Asymptotic confidence intervale for the median in the
                    ## one-sample case.  When used with paired values this
                    ## gives a confidence interval for mean(x) - mean(y).
                    ## Algorithm not published, thus better documented here.
                    x <- x + mu
                    alpha <- 1 - conf.level
                    ## These are sample based limits for the median
                    mumin <- min(x)
                    mumax <- max(x)
                    ## wdiff(d, zq) returns the absolute difference between  
                    ## the asymptotic Wilcoxon statistic of x - mu - d and
                    ## the quantile zq.
                    CORRECTION.CI <- 0
                    wdiff <- function(d, zq) {
                        xd <- x - d
                        xd <- xd[xd != 0]
                        nx <- length(xd)
                        dr <- rank(abs(xd))
                        zd <- sum(dr[xd > 0])
                        zd <- (zd - wmean) / sqrt(wvar)
                        zd - zq
                    }
                    ## Here we search the root of the function wdiff on the set
                    ## c(mumin, mumax).
                    ##
                    ## This returns a value from c(mumin, mumax) for which
                    ## the asymptotic Wilcoxon statistic is equal to the
                    ## quantile zq.  This means that the statistic is not
                    ## within the critical region, and that implies that d
                    ## is a confidence limit for the median.
                    ##
                    ## As in the exact case, interchange quantiles.
                    cint <- switch(alternative, "two.sided" = {
                        l <- uniroot(wdiff, c(mumin, mumax), tol=1e-4,
                                     zq=qnorm(alpha/2, lower.tail=FALSE))$root
                        u <- uniroot(wdiff, c(mumin, mumax), tol=1e-4,
                                     zq=qnorm(alpha/2))$root
                        c(l, u)
                    }, "greater"= {
                        l <- uniroot(wdiff, c(mumin, mumax), tol=1e-4,
                                     zq=qnorm(alpha, lower.tail=FALSE))$root
                        c(l, +Inf)
                    }, "less"= {
                        u <- uniroot(wdiff, c(mumin, mumax), tol=1e-4,
                                     zq=qnorm(alpha))$root
                        c(-Inf, u)
                    })
                    attr(cint, "conf.level") <- conf.level    
                    ESTIMATE <- uniroot(wdiff, c(mumin, mumax), tol=1e-4,
                                    zq=0)$root
                    names(ESTIMATE) <- "(pseudo)median"
                 }

	}
    }
    else {
        if(length(y) < 1) 
            stop("not enough y observations")
        METHOD <- "Exact Wilcoxon rank sum test"
        r <- rank(c(x - mu, y))
        n.x <- length(x)
        n.y <- length(y)
        if(is.null(exact)) 
            exact <- (n.x < 50) && (n.y < 50)
        STATISTIC <- sum(r[seq(along = x)]) - n.x * (n.x + 1) / 2
        names(STATISTIC) <- "W"
        TIES <- (length(r) != length(unique(r)))
        if(exact) {
            PVAL <-
                switch(alternative,
                       "two.sided" = pperm(STATISTIC + n.x*(n.x +1)/2, r,
                                           n.x, alternative="two.sided", pprob=TRUE),
                       "greater" = pperm(STATISTIC + n.x*(n.x+1)/2, r,
                                       n.x, alternative="greater", pprob=TRUE), 
                           "less" = pperm(STATISTIC+ n.x*(n.x +1)/2, r, n.x,
                                          alternative="less", pprob=TRUE))
             MIDP <- PVAL$PPROB    
             PVAL <- PVAL$PVALUE
             if(conf.int) {
                 ## Exact confidence interval for the location parameter 
                 ## mean(y) - mean(x) in the two-sample case (cf. the
                 ## one-sample case).
                 if (mu != 0) r <- rank(c(x,y))
                 alpha <- 1 - conf.level
                 diffs <- sort(outer(x, y, "-"))
                 if (TIES) {
                   fs <- function(d)
                     sum(rank(c(x-d,y))[seq(along = x)]) - n.x * (n.x + 1) / 2
                   w <- sapply(diffs, fs)
                 } else {
                   w <- (n.x*n.y):1
                 }
                 cint <-
                     switch(alternative,
                            "two.sided" = {
                                qu <- qperm(alpha/2, r, n.x) - n.x*(n.x+1)/2
                                ql <- qperm(1-alpha/2, r, n.x) - n.x*(n.x+1)/2
                                if (qu <= min(w)) lci <- max(diffs)
                                else lci <- min(diffs[w <= qu])
                                if (ql >= max(w)) uci <- min(diffs)
                                else uci <- max(diffs[w > ql])
                                c(uci, lci)
                            },
                            "greater"= {
                                ql <- qperm(1-alpha, r, n.x) - n.x*(n.x+1)/2
                                if (ql >= max(w)) uci <- min(diffs)
                                else uci <- max(diffs[w > ql])
                                c(uci, +Inf)
                             },
                             "less"= {
                                 qu <- qperm(alpha, r, n.x) - n.x*(n.x+1)/2
                                 if (qu <= min(w)) lci <- max(diffs)
                                 else lci <- min(diffs[w <= qu])
                                 c(-Inf, lci)
                              })
                  attr(cint, "conf.level") <- conf.level    
                  wmean <- n.x/(n.x+n.y)*sum(r) - n.x*(n.x+1)/2
                  ESTIMATE <- mean(c(min(diffs[w <= ceiling(wmean)]),
                                     max(diffs[w > wmean])))
                  names(ESTIMATE) <- "difference in location"
             }
        } else {
            METHOD <- "Asymptotic Wilcoxon rank sum test"
            N <- n.x + n.y
            wmean <- n.x/N*sum(r)
            wvar <- n.x*n.y/(N*(N-1))*sum((r - wmean/n.x)^2)
            PVAL <- pnorm((STATISTIC + n.x*(n.x+1)/2 - wmean)/sqrt(wvar))
            if (alternative == "greater")
               PVAL <- 1 - PVAL
            if(alternative == "two.sided") 
                PVAL <- 2 * min(PVAL, 1 - PVAL)
            if(conf.int) {
                ## Asymptotic confidence interval for the location
                ## parameter mean(x) - mean(y) in the two-sample case
                ## (cf. one-sample case).
                ##
                ## Algorithm not published, for a documentation see the
                ## one-sample case.
                alpha <- 1 - conf.level
                mumin <- min(x) - max(y)
                mumax <- max(x) - min(y)
                CORRECTION.CI <- 0
                wdiff <- function(d, zq) {
                    dr <- rank(c(x - d, y))
                    dz <- sum(dr[seq(along = x)])
                    dz <- (dz - wmean) / sqrt(wvar)
                    dz - zq
                }
                cint <- switch(alternative, "two.sided" = {
                    l <- uniroot(wdiff, c(mumin, mumax), tol=1e-4,
                                  zq=qnorm(alpha/2, lower.tail=FALSE))$root
                    u <- uniroot(wdiff, c(mumin, mumax), tol=1e-4,
                                  zq=qnorm(alpha/2))$root
                    c(l, u)
                }, "greater"= {
                    l <- uniroot(wdiff, c(mumin, mumax), tol=1e-4,
                                  zq=qnorm(alpha, lower.tail=FALSE))$root
                    c(l, +Inf)
                }, "less"= {
                    u <- uniroot(wdiff, c(mumin, mumax), tol=1e-4,
                                  zq=qnorm(alpha))$root
                    c(-Inf, u)
                })
                attr(cint, "conf.level") <- conf.level    
                ESTIMATE <- uniroot(wdiff, c(mumin, mumax), tol=1e-4,
                                    zq=0)$root
                names(ESTIMATE) <- "difference in location"
            }
        }
    }

    if (!is.null(MIDP)) {
        names(MIDP) <- "point prob"
        RVAL <- list(statistic = STATISTIC,
                     pointprob = MIDP,
                     p.value = PVAL, 
                     null.value = c(mu = mu),
                     alternative = alternative,
                     method = METHOD, 
                     data.name = DNAME)
    } else {
        RVAL <- list(statistic = STATISTIC,
                     p.value = PVAL, 
                     null.value = c(mu = mu),  
                     alternative = alternative,
                     method = METHOD, 
                     data.name = DNAME)
    }
    if(conf.int) {
        RVAL$conf.int <- cint
        RVAL$estimate <- ESTIMATE
    }
    class(RVAL) <- "htest"
    return(RVAL)
}

wilcox.exact.formula <-
function(formula, data, subset, na.action, ...)
{
    if(missing(formula)
       || (length(formula) != 3)  
       || (length(attr(terms(formula[-2]), "term.labels")) != 1)
       || (length(attr(terms(formula[-3]), "term.labels")) != 1))
        stop("formula missing or incorrect")
    if(missing(na.action))
        na.action <- getOption("na.action")
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
    DATA <- split(mf[[response]], g)
    names(DATA) <- c("x", "y")
    y <- do.call("wilcox.exact", c(DATA, list(...)))
    y$data.name <- DNAME
    y
}
