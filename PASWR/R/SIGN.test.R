SIGN.test <-
function(x, y = NULL, md = 0, alternative = "two.sided", conf.level = 0.95)
{
    choices <- c("two.sided", "greater", "less")
    alt <- pmatch(alternative, choices)
    alternative <- choices[alt]
    if(length(alternative) > 1 || is.na(alternative))
        stop("alternative must be one \"greater\", \"less\", \"two.sided\"")
    if(!missing(md))
        if(length(md) != 1 || is.na(md))
            stop("median must be a single number")
    if(!missing(conf.level))
        if(length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
            stop("conf.level must be a number between 0 and 1")

    if( is.null(y) )
    {
    # One-Sample Sign-Test Exact Test
    dname <- paste(deparse(substitute(x)))
    x <- sort(x)
    diff <- (x - md)
    n <- length(x)
    nt <- length(x) - sum(diff == 0)
    s <- sum(diff > 0)
    estimate <- median(x)
    method <- c("One-sample Sign-Test")
    names(estimate) <- c("median of x")
    names(md) <- "median"
    names(s) <- "s"
    CIS <- "Conf Intervals"
        if(alternative == "less")
        {
        # zobs <- (s-0.5*n)/sqrt(n*0.25)
        pval <- sum(dbinom(0:s, nt, 0.5))
        # Note: Code uses linear interpolation to arrive at the confidence intervals.
        loc <- c(0:n)
        prov <- (dbinom(loc, n, 0.5))
        k <- loc[cumsum(prov) > (1 - conf.level)][1]

            if(k < 1)
            {
            conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
            xl <- -Inf
            xu <- x[n]
            ici <- c(xl, xu)
            }
            else
            {
            ci1 <- c(-Inf, x[n - k + 1])
            acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
            ci2 <- c(-Inf, x[n - k])
            acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
            xl <- -Inf
            xu <- (((x[n - k + 1] - x[n - k]) * (conf.level - acl2))/(acl1 - acl2)) + x[n - k]
            ici <- c(xl, xu)
            }

        }
        else if(alternative == "greater")
        {
        pval <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
        loc <- c(0:n)
        prov <- (dbinom(loc, n, 0.5))
        k <- loc[cumsum(prov) > (1 - conf.level)][1]

            if(k < 1)
            {
            conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
            xl <- x[1]
            xu <- Inf
            ici <- c(xl, xu)
            }
            else
            {
            ci1 <- c(x[k], Inf)
            acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
            ci2 <- c(x[k + 1], Inf)
            acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
            xl <- (((x[k] - x[k + 1]) * (conf.level - acl2))/(acl1 - acl2)) + x[k + 1]
            xu <- Inf
            ici <- c(xl, xu)
            }
        }
        else
        {
        p1 <- sum(dbinom(0:s, nt, 0.5))
        p2 <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
        pval <- min(2 * p1, 2 * p2, 1)
        loc <- c(0:n)
        prov <- (dbinom(loc, n, 0.5))
        k <- loc[cumsum(prov) > (1 - conf.level)/2][1]

            if(k < 1)
            {
            conf.level <- (1 - 2 * (sum(dbinom(k, n, 0.5))))
            xl <- x[1]
            xu <- x[n]
            ici <- c(xl, xu)
            }
            else
            {
            ci1 <- c(x[k], x[n - k + 1])
            acl1 <- (1 - 2 * (sum(dbinom(0:k - 1, n, 0.5))))
            ci2 <- c(x[k + 1], x[n - k])
            acl2 <- (1 - 2 * (sum(dbinom(0:k, n, 0.5))))
            xl <- (((x[k] - x[k + 1]) * (conf.level - acl2))/(acl1 - acl2)) + x[k + 1]
            xu <- (((x[n - k + 1] - x[n - k]) * (conf.level - acl2))/(acl1 - acl2)) + x[n - k]
            ici <- c(xl, xu)
            }

        }
    }



    else
    {
    #   Paired-Samples Sign Test
    if(length(x)!=length(y))
        stop("Length of x must equal length of y")
    xy <- sort(x-y)
    diff <- (xy - md)
    n <- length(xy)
    nt <- length(xy) - sum(diff == 0)
    s <- sum(diff > 0)
    dname <-  paste(deparse(substitute(x)), " and ", deparse(substitute(y)), sep = "")
    estimate <- median(xy)
    method <- c("Dependent-samples Sign-Test")
    names(estimate) <- c("median of x-y")
    names(md) <- "median difference"
    names(s) <- "S"
    CIS <- "Conf Intervals"
    if(alternative == "less")
    {
        pval <- sum(dbinom(0:s, nt, 0.5))
        # Note: Code uses linear interpolation to arrive at the confidence intervals.
        loc <- c(0:n)
        prov <- (dbinom(loc, n, 0.5))
        k <- loc[cumsum(prov) > (1 - conf.level)][1]

            if(k < 1)
            {
            conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
            xl <- -Inf
            xu <- xy[n]
            ici <- c(xl, xu)
            }
            else
            {
            ci1 <- c(-Inf, xy[n - k + 1])
            acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
            ci2 <- c(-Inf, xy[n - k])
            acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
            xl <- -Inf
            xu <- (((xy[n - k + 1] - xy[n - k]) * (conf.level - acl2))/(acl1 - acl2)) + xy[n - k]
            ici <- c(xl, xu)
            }

    }
    else if(alternative == "greater")
    {
        pval <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
        loc <- c(0:n)
        prov <- (dbinom(loc, n, 0.5))
        k <- loc[cumsum(prov) > (1 - conf.level)][1]

        if(k < 1)
            {
            conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
            xl <- xy[1]
            xu <- Inf
            ici <- c(xl, xu)
            }
            else
            {
            ci1 <- c(xy[k], Inf)
            acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
            ci2 <- c(xy[k + 1], Inf)
            acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
            xl <- (((xy[k] - xy[k + 1]) * (conf.level - acl2))/(acl1 - acl2)) + xy[k + 1]
            xu <- Inf
            ici <- c(xl, xu)
            }

     }
        else
        {
        p1 <- sum(dbinom(0:s, nt, 0.5))
        p2 <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
        pval <- min(2 * p1, 2 * p2, 1)
        loc <- c(0:n)
        prov <- (dbinom(loc, n, 0.5))
        k <- loc[cumsum(prov) > (1 - conf.level)/2][1]
        if(k < 1)
            {
            conf.level <- (1 - 2 * (sum(dbinom(k, n, 0.5))))
            xl <- xy[1]
            xu <- xy[n]
            ici <- c(xl, xu)
            }
            else
            {
            ci1 <- c(xy[k], xy[n - k + 1])
            acl1 <- (1 - 2 * (sum(dbinom(0:k - 1, n, 0.5))))
            ci2 <- c(xy[k + 1], xy[n - k])
            acl2 <- (1 - 2 * (sum(dbinom(0:k, n, 0.5))))
            xl <- (((xy[k] - xy[k + 1]) * (conf.level - acl2))/(acl1 - acl2)) + xy[k + 1]
            xu <- (((xy[n - k + 1] - xy[n - k]) * (conf.level - acl2))/(acl1 - acl2)) + xy[n - k]
            ici <- c(xl, xu)
            }
        }
    }

    # Below is how I had the stuff print at first...showing the exact intervals as
    # well as the interpolated interval.  Changed format to standard htest class.
    # Test.Values <- t(as.matrix(c(s, pval)))
    # values <- c("S", "p-value")
    # dimnames(Test.Values) <- list(NULL, values)

    if(k < 1)
        {
        cint <- ici
        attr(cint, "conf.level") <- conf.level
        rval <- structure(list(statistic = s, p.value = pval, estimate = estimate, null.value = md,
        alternative = alternative, method = method, data.name = dname, conf.int=cint ))
        oldClass(rval) <- "htest"
        return(rval)
        }
    else
    {
    result1 <- c(acl2, ci2)
     result2 <- c(conf.level, ici)
     result3 <- c(acl1, ci1)
     Confidence.Intervals <- round(as.matrix(rbind(result1, result2, result3)),4)
     cnames <- c("Conf.Level", "L.E.pt", "U.E.pt")
     rnames <- c("Lower Achieved CI", "Interpolated CI", "Upper Achieved CI")
     dimnames(Confidence.Intervals) <- list(rnames, cnames)
    # return(Test.Values, Confidence.Intervals)
    cint <- ici

    attr(cint, "conf.level") <- conf.level

    rval <- structure(list(statistic = s, parameter=NULL, p.value = pval,
      conf.int = cint, estimate = estimate, null.value = md,
      alternative = alternative, method = method, data.name = dname ))
    ## Returns both htest format and matrix of CIs
    ## attr(rval, "class") <- "htest"
        oldClass(rval) <- "htest"
    ##   return(rval, Confidence.Intervals)
    # return(rval)
     print(rval)
     return(Confidence.Intervals)
    }
}

