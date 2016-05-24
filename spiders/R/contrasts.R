##' @title linear contrast of c_st
##'
##' @description estimates linear contrasts of the elements of c, c_s, c_t, or c_st from a \code{\link{predPref}} object
##'
##' @return A list with class '"htest"' containing the following components:
##' 
##' statistic: the value of the t-statistic.
##' 
##' parameter: the degrees of freedom for the t-statistic.
##' 
##' p.value: the p-value for the test.
##' 
##' conf.int: a confidence interval for the mean appropriate to the specified
##' alternative hypothesis.
##' 
##' estimate: the estimated mean or difference in means depending on whether
##' it was a one-sample test or a two-sample test.
##' 
##' null.value: the specified hypothesized value of the mean or mean difference
##' depending on whether it was a one-sample test or a two-sample test.
##' 
##' alternative: a character string describing the alternative hypothesis.
##' 
##' method: a character string indicating what type of t-test was performed.
##' 
##' data.name: a character string giving the names of the data.
##' 
##' @details The input vector b performs the linear transformation
##' t(b) \%*\% matrix(c_st), so that c_st becomes a column vector by indexing
##' t first and then s.  Hence there is no requirement of a linear
##' contrast, any linear transformation such that
##' t(b) \%*\% matrix(1, nrow=length(b)) != 0 is allowed.
##'
##' Of the two estimated hypotheses in the underlying call
##' to \code{\link{predPref}}, the linear transformation b is applied to the
##' hypothesis that is determined by the choice of \code{sig.level}.
##'
##' @param x a predPref object as fit by the eponymous function
##' @param b a vector to linearly transform c_st
##' @param mu a number to test the linear contrast against in the null
##' @param alternative string to specify alternative hypothesis
##' @param conf.level confidence level of the interval
##' @param sig.level determines null/alternative hypothesis value of c_st from predPref
##' 
##' @examples
##' # set parameters
##' Predators <- Traps <- 100
##' PreySpecies <- 2
##' Times <- 5
##' g <- matrix(sqrt(2), nrow=Times, ncol=PreySpecies)     # gamma
##' l <- matrix(seq(0.4,1.8,length.out=5)*sqrt(2), nrow=Times, ncol=PreySpecies) # ct
##'
##' # fit model and contrast
##' \dontrun{
##' set.seed(0)
##' fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=FALSE) # p-value=0.305
##' pref <- predPref(fdata$eaten, fdata$caught, hypotheses=c('ct', 'cst'))
##' testC(pref, b = c(0,1, -1, 0, 0)) # p-value > sig.level => ct is used, not cst
##' }
##' 
##' @export
testC <- function(x, b, mu = 0, alternative = c("two.sided", "less", "greater"), conf.level = 0.95, sig.level=0.05) {
    
    ## some check on input; stolen from t.test
    if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
        stop("'mu' must be a single number")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                                 conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")

    ## get appropriate estimates
    if ( x$LRT$p.value < sig.level ) {
            C <- matrix(x$alt$c)
            lenC <- nrow(C); lc <- seq_len(lenC)
            varC <- x$alt$var[lc,lc]
    } else {
        if ( !is.null(x$null$c) ) {
            C <- matrix(x$null$c)
            lenC <- nrow(C); lc <- seq_len(lenC)
            varC <- x$null$var[lc,lc]
        } else {
            stop("Don't yet know how to handle 1 = lambda_{st}/gamma_{st}.")
        }        
    }

    ## checks on contrast
    if ( length(b) != lenC ) {
        stop("length of contrast b and length of C don't match.")
    }
        

    ## calculations
    stat <- t(b)%*%C
    scale <- t(b)%*%varC%*%b
    stdev <- sqrt(scale)

    ## test and confidence interval
    alpha <- 1-conf.level
    a2 <- alpha/2
    alternative <- match.arg(alternative)
    if ( alternative == "two.sided" ) {
        pval <- 2*pnorm(-abs(stat), mean=mu, sd=stdev)
        interval <- qnorm(c(a2, 1-a2), mean=stat, sd=stdev)
    } else if ( alternative == "less" ) {
        pval <- pnorm(stat, mean=mu, sd=stdev)
        interval <- c(-Inf, qnorm(conf.level, mean=stat, sd=stdev))
    } else if ( alternative == "greater" ) {
        pval <- pnorm(stat, mean=mu, sd=stdev, lower.tail=F)
        interval <- c(qnorm(alpha, mean=stat, sd=stdev), Inf)
    }

    ## info to fill in print method of class htest
    names(stat) <- "b^t * C"
    names(mu) <- "linear contrast"
    method <- paste("Linear Contrast: ", paste(t(b), collapse=" "), sep="")
    attr(interval, "conf.level") <- conf.level
    est <- as.vector(C)
    names(est) <- paste("mean of c_", lc, sep="")

    out <- list(statistic = stat, p.value = pval,
                conf.int = interval, null.value = mu,
                alternative = alternative, method = method,
                estimate = est, data.name=x$data.name)
    class(out) <- "htest"
    out
}
