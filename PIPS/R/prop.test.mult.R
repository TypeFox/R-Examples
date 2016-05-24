###########################################################################
## PROGRAM:   prop.test.mult.R
## PURPOSE:
##
## Function to perform a 1 or 2 sample test of proportions.  The code
## comes largely from prop.test in the stats package, but differs in
## that it allows multiple unrelated tests to be carried out
## simultaneously.  If x2 is not NULL (i.e., we're in the 2-sample
## case), then x1 and x2 should be vectors of the same length.  The
## number of tests performed is always the same as the length of x1.
## The arguments n1 and n2 should either be the same length as x1 and
## x2, respectively, or they should each have a single value which
## will be used to for all tests.
## 
## INPUT PARAMETERS: 
## x1, x2: Number of successes from group 1,2
## n1, n2: Number of group 1, group 2 trials.  This is either a 
##  vector of the same length as x1 or a single number used as the 
##  group 1 size for all the tests.
## samples: Can pass samples==1 to force a one-sample test.  A 
##  one-sample test is also performed if x2 is null. 
## alternative: Alternative hypothesis
## conf.level: Confidence level (0-1). Default is 0.95. 
## correct: TRUE/FALSE - use continuity correction?  Default is TRUE.
## as.vector: TRUE/FALSE - If TRUE and only one test was done 
##  (length(x1)==1) then return results as a vector (this is more 
##  convenient for the calling program.  Otherwise, results will be
##  returned as a matrix
## 
## OUTPUT (RETURN VALUE):
## The function returns either a vector (if as.vector=TRUE - with elements
## {estimate,upper, lower}) or a matrix (with columns {estimate,upper,lower} and
## one row for every test.
## MACROS USED:   None
## CALLED BY:     pred.int.R
## AUTHOR:        Daniel Muenz 
## CREATION DATE: 2010
## NOTES:
## MODIFICATIONS:
## [RG20120107] Ray Griner standardized stop messages
## [RG20120830] Ray Griner stdize header
######################################################################
prop.test.mult <-
function(x1, n1, x2=NULL, n2=NULL, samples=NULL, p=NULL,
         alternative=c("two.sided", "less", "greater"),
         conf.level=0.95, correct=TRUE, as.vector=FALSE)
{
    if (length(n1)!=1 && length(n1)!=length(x1)) 
      stop ("'n1' must have the same length as 'x1' or have length 1") ## [RG20120107]

    if (!is.null(x2) && length(n2)!=1 && length(n2)!=length(x2)) 
      stop ("if 'x2' is not null, 'n2' must have the same length as 'x2' or have length 1") ## [RG20120107] 

    if (samples == 1 || is.null(x2)) k <- 1
    else k <- 2

    if (is.null(p) && (k == 1)) p <- 0.5

    alternative <- match.arg(alternative)

    if ((k == 2) && !is.null(p))
        alternative <- "two.sided"

    correct <- as.logical(correct)

    ######################################################################
    ## Daniel wrote code implementing the Yates continuity correction but
    ##  pred.int always calls this with correct==FALSE, so Im not testing
    ##  it now. If someone needs this correction, I will test this and
    ##  remove the error message
    ######################################################################
    if (correct==TRUE) stop('yates continuity correction not supported') ## [RG20120107] 

    NVAL <- p
    CINT <- NULL
    YATES <- ifelse(correct, 0.5, 0)
    z <- ifelse(alternative == "two.sided", qnorm((1 + conf.level)/2),
                qnorm(conf.level))

    if (k == 1) {
        ESTIMATE <- x1 / n1
        YATES <- sapply(abs(x1 - n1*p), min, YATES)
        z22n <- z^2 / (2*n1)

        ## Upper bound
        p.c <- ESTIMATE + YATES/n1
        ######################################################################
        # These next two lines are slightly different from the R fn prop.test??
        ######################################################################
        p.c <- sapply(p.c, min, 1)
        p.u <- (p.c + z22n + z * sqrt(p.c*(1 - p.c)/n1 + z22n/(2*n1))) / (1 + 2*z22n)
        p.u <- sapply(p.u, min, 1)

        ## Lower bound
        p.c <- ESTIMATE - YATES/n1
        ######################################################################
        # These next two lines are slightly different from the R fn prop.test??
        ######################################################################
        p.c <- sapply(p.c, max, 0)
        p.l <- (p.c + z22n - z * sqrt(p.c*(1 - p.c)/n1 + z22n/(2*n1))) / (1 + 2*z22n)
        p.l <- sapply(p.l, max, 0)

        CINT <- switch(alternative,
                       two.sided = cbind(p.l, p.u),
                       greater = cbind(p.l, 1),
                       less = cbind(0, p.u))
    }
    else if ((k == 2) & is.null(p)) {
        p1 <- x1 / n1
        p2 <- x2 / n2
        ESTIMATE <- p1 - p2
        YATES <- sapply((abs(ESTIMATE)/(1/n1 + 1/n2)), min, YATES)

        WIDTH <- z * sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2) + YATES*(1/n1 + 1/n2)
        LOWER <- if (alternative == "less") -1
                 else sapply(ESTIMATE - WIDTH, max, -1)
        UPPER <- if (alternative == "greater") 1
                 else sapply(ESTIMATE + WIDTH, min, 1)
        CINT <- cbind(LOWER, UPPER)
    }

    rval <- cbind(ESTIMATE, CINT)
    if (as.vector && length(x1)==1) {
        rval <- as.vector(rval);
        names(rval) <- c("point","lower","upper")
    }
    else
        dimnames(rval)[[2]] <- c("point","lower","upper")
    return(rval)
}
