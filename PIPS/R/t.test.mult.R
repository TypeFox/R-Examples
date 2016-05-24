###########################################################################
## PROGRAM: t.test.mult.R
## PURPOSE:
##
## Function to perform a 1 or 2 sample t-test.  The code comes largely
##  from t.test in the stats package, but differs in that it is used
##  with summary statistics rather than raw data, and it allows
##  multiple unrelated tests to be carried out simultaneously.
##
## Note that several of the input parameters (mean1, var1, mean2, var2) can 
##  accept vectors as input.  All these vectors must have the same length.
##  If these are vectors, then a t-test will be conducted for each element
##  of the vectors.  The other parameters will never be vectors.  For example
##  n1, n2 always represent the sample size of the two groups, and the 
##  program uses this same sample size for every t-test it performs.
## 
## INPUT PARAMETERS:
## mean1,mean2 - Mean for Group 1/2.  Either a single value or a vector.
## var1,var2 - Variance for Group 1/2.  Either a single value or a vector.
## n1,n2 - Number in Group 1/2.  Each of these is a single value, not a 
##  vector. It is assumed all t-tests have the same numbers in the groups
## samples - Number of groups we are comparing.  If 1, a one-sample t-test
##  is performed.  Otherwise, a 2 sample t-test is done.
## var.equal - TRUE/FALSE whether to assume variances in the two groups are 
##  equal
## conf.level - Confidence level (between 0 and 1).  Default is 0.95
## as.vector - TRUE/FALSE. Default is FALSE.  If TRUE and length(mean1)=1, the
##  result is output as a vector.  The function is called using TRUE when
##  doing a t-test on the observed data for each comparison, because it's 
##  more convenient for the calling function to have the output data as a vector
## 
## OUTPUT (RETURN VALUE):
## The function returns either a vector (if as.vector=TRUE - with elements
## {estimate,upper, lower}) or a matrix (with columns {estimate,upper,lower} and
## one row for every test.
##
## MACROS USED: None
## CALLED BY:   None
## AUTHOR: Daniel Muenz
## CREATION DATE: 2010
## NOTES:
## MODIFICATIONS:
## [RG20120107] Ray Griner standardized stop messages
## [RG20120830] Ray Griner standardized program header
###########################################################################
t.test.mult <-
function(mean1, var1, n1, mean2=NULL, var2=NULL, n2=NULL,
         samples=NULL, alternative=c("two.sided", "less", "greater"),
         mu=0, var.equal=FALSE, conf.level=0.95, as.vector=FALSE)
{
    alternative <- match.arg(alternative)

    #################################################################
    ## Check that mean1, var1, mean2, var2 all have the same length #
    #################################################################
    ## [RG20120107]S
    if (length(mean1) != length(var1)) {
      stop("'mean1' and 'var1' must have the same length.")
    }
    if (any(is.na(var1))) { stop("'var1' must not be null or NA") }

    if (!is.null(mean2) && (length(mean2) != length(var2))) {
      stop("'mean2' and 'var2' must have the same length.")
    }
    if (!is.null(mean2) && any(is.na(var2))) { stop("'var2' must not be null or NA if 'mean2' is not null") }

    if (!is.null(mean2) && length(mean1)!=length(mean2)) {
      stop("'mean1' and 'mean2' must have the same length if both are not null")
    }
    if (length(n1)!=1) { stop("'n1' must have length 1") }
    if ((length(n2)!=1) && (length(n2)!=NA)) { stop("'n2' must be null or have length 1") }
    ## [RG20120107]E

    if (samples == 1) {
        estimate <- mean1
	df <- rep(n1-1, length(mean1))
	stderr <- sqrt(var1/n1)
	tstat <- (mean1-mu)/stderr
    }
    else {
	estimate <- mean1 - mean2
	if (var.equal) {
            #######################################
            # This looks OK. - Ray                #
            #######################################
	    df <- rep(n1+n2-2, length(mean1))
            v <- 0
            if(n1 > 1) v <- v + (n1-1)*var1
            if(n2 > 1) v <- v + (n2-1)*var2
	    v <- v/df
	    stderr <- sqrt(v*(1/n1+1/n2))
	}
        else {
            ######################################################
            # Using Sattherthwaites approximation. Also looks OK #
            ######################################################
	    stderr1 <- sqrt(var1/n1)
	    stderr2 <- sqrt(var2/n2)
	    stderr <- sqrt(stderr1^2 + stderr2^2)
	    df <- stderr^4/(stderr1^4/(n1-1) + stderr2^4/(n2-1))
	}
        tstat <- (estimate - mu)/stderr
    }

    #####################################################
    # Get the confidence intervals 
    #####################################################
    if (alternative == "less") {
	cint <- cbind(-Inf, tstat + qt(conf.level, df) )
    }
    else if (alternative == "greater") {
	cint <- cbind(tstat - qt(conf.level, df), Inf)
    }
    else {
        cint <- qt((1+conf.level)/2, df)
	cint <- tstat + cbind(-cint, cint)
    }
    cint <- mu + cint * stderr

    rval <- cbind(estimate, cint)

    ######################################################
    # If we only calculated one t-test, output it as     #  
    #  a vector.  Otherwise, output it as a matrix where #
    #  the columns are the estimate, lower, and upper    #
    #  bounds, and the rows are each t-test that we did  #
    ######################################################
    if (as.vector && length(mean1)==1) {
        rval <- as.vector(rval);
        names(rval) <- c("point","lower","upper")
    }
    else
        dimnames(rval)[[2]] <- c("point","lower","upper")
    return(rval)
}
