#############################################################################
##                                                                         ##
##   Tests for Runuran functions                                           ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Remark: You must use named arguments when calling the test routines!  ##
##                                                                         ##
#############################################################################

## --- Run tests? -----------------------------------------------------------

unur.run.tests <- TRUE   ## <-- change to FALSE if tests should not be performed

if (!unur.run.tests) {
        cat("\nRunuran tests not performed!\n\n")
        quit(save="no",status=0,runLast=FALSE)
}

## --- Test Parameters ------------------------------------------------------

## size of sample for test
samplesize <- 1.e5

## level of significance
alpha <- 1.e-3

## number of repetitions
n.rep.domains <- 2
n.rep.params  <- 5

## --- Global counters ------------------------------------------------------

unur.envir <- new.env()
assign("pvals", numeric(0), env = unur.envir)
assign("n.warns", 0, env = unur.envir)

## --- Load library ---------------------------------------------------------

library(Runuran)

#############################################################################
##                                                                          #
##  Error Measures for Inversion Algorithms                                 #
##                                                                          #
#############################################################################

## --- U-error --------------------------------------------------------------

unur.cont.uerror <- function (n, aqfunc, pfunc) {
  ## Compute maximal u-error(u) = | u - pfunc( aqfunc (u) ) |
  ##
  ##  n      ... sample size
  ##  aqfunc ... approximate quantile function
  ##  pfunc  ... cumulative probability function (CDF)
  ##

  ## check input
  if (missing(aqfunc) || missing(pfunc))
    stop ("aqfunc and pfunc required!")

  ## u-values (equidistributed)
  u <- (0:(n-1))/n + 1/(2*n)
  
  ## u-error
  ue <- abs( u - pfunc( aqfunc (u) ) )

  ## return maximal u-error
  max(ue)
}

## --- X-error --------------------------------------------------------------

unur.xerror <- function (n, aqfunc, qfunc) {
  ## Compute maximal x-error(u) = | aqfunct(u) - qfunc(u) |
  ##
  ##  n      ... sample size
  ##  aqfunc ... approximate quantile function
  ##  qfunc  ... (exact) quantile function
  ##

  ## check input
  if (missing(aqfunc) || missing(qfunc))
    stop ("aqfunc and qfunc required!")

  ## u-values (equidistributed)
  u <- (0:(n-1))/n + 1/(2*n)
  
  ## u-error
  xe <- abs( qfunc(u) - aqfunc(u) )

  ## return maximal x-error
  max(xe)
}


#############################################################################
##                                                                          #
##  Continuous univariate Distributions                                     #
##                                                                          #
#############################################################################

## --- CONT: Function for running chi^2 goodness-of-fit test ----------------

unur.test.cont <- function (distr, rfunc=NULL, pfunc=NULL, domain, ...) {
        ##  Run a chi^2 test and evaluate p-value.
        ##  Repeat test once if it fails
        ##  (we do not check the validity of the algorithm
        ##   but the correct implementation.)
        ##
        ##  distr  ... name of distribution (as used for [rpqd]... calls
        ##  rfunc  ... random number generation function
        ##  pfunc  ... probability function (CDF)
        ##  domain ... domain of distribution
        ##  '...'  ... list of parameters of distribution
        ##

        ## -- domain ?
        have.domain = ifelse( missing(domain), FALSE, TRUE )
        lb <- ifelse( have.domain, domain[1], -Inf)
        ub <- ifelse( have.domain, domain[2], Inf)
        
        ## -- text for calling function
        cat(distr,"(",paste(...,sep=","),")",
            ifelse( have.domain,paste(" domain=(",signif(lb),",",signif(ub),"): ",sep=""),": "),
            sep="")

        ## -- function for generating a random sample
        if (is.null(rfunc)) {
                rfuncname <- paste("ur",distr,sep="")
                if (!exists(rfuncname))
                        stop("undefined function '",rfuncname,"'")
                rfunc <- match.fun(rfuncname, descend=FALSE)
        }
        ## -- function for computing CDF
        if (is.null(pfunc)) {
                pfuncname <- paste("p",distr,sep="")
                if (!exists(pfuncname))
                        stop("undefined function '",pfuncname,"'")
                pfunc <- match.fun(pfuncname, descend=FALSE)
        }
        
        ## -- run test and repeat once when failed the first time
        for (i in 1:2) {

                ## -- random sample
                if (have.domain)
                        x <- rfunc(samplesize,lb=lb,ub=ub,...)
                else
                        x <- rfunc(samplesize,...)
                ## -- test domain (if given)
                if (have.domain) {
                        too.small <- length(x[x<lb])
                        too.large <- length(x[x>ub])
                        if (too.small > 1 || too.small > 1) {
                                too.small <- 100*too.small/samplesize
                                too.large <- 100*too.large/samplesize
                                stop("X out of domain (",too.small,"%|",too.large,"%)", call.=FALSE)
                        }
                }
                
                ## -- transform into uniform distribution
                u <- pfunc(x,...)
                ## -- make histogram of with equalsized bins (classified data)
                nbins <- as.integer(sqrt(samplesize))
                breaks <- (0:nbins)/nbins
                if (have.domain) {
                        u.lb = pfunc(lb,...)
                        u.ub = pfunc(ub,...)
                        breaks <- u.lb + breaks*(u.ub-u.lb)
                        breaks[length(breaks)] <- u.ub
                }
                h <- hist(u,plot=F,breaks=breaks)$count
                ## -- run unur.chiq.test --
                pval <- chisq.test(h)$p.value
                ## -- check p-value
                if (pval > alpha) { # test passed
                        message("chisq test PASSed with p-value=",signif(pval))
                        break
                }

                ## -- test failed
                if (i>1) { # second try failed
                        stop("chisq test FAILED!  p-value=",signif(pval), call.=FALSE)
                }
                else {
                        warning("first chisq test FAILed with p-value=",signif(pval),
                                call.=FALSE,immediate.=TRUE)
                        assign("n.warns", get("n.warns",envir=unur.envir)+1, env=unur.envir)

                }
        }

        ## -- update list of p-values
        assign("pvals", append(get("pvals",envir=unur.envir), pval), env=unur.envir)

} ## --- end of unur.test.cont() ---


#############################################################################
##                                                                          #
##  Discrete univariate Distributions                                       #
##                                                                          #
#############################################################################

## --- DISCR: Function for running chi^2 goodness-of-fit test ---------------

unur.test.discr <- function (distr, rfunc=NULL, dfunc=NULL, pv, domain, ...) {
        ##  Run a chi^2 test and evaluate p-value.
        ##  Repeat test once if it fails
        ##  (we do not check the validity of the algorithm
        ##   but the correct implementation.)
        ##
        ##  distr  ... name of distribution (as used for [rpqd]... calls
        ##  rfunc  ... random number generation function
        ##  dfunc  ... probability mass function
        ##  pv     ... probability vector
        ##  domain ... domain of distribution
        ##  '...'  ... list of parameters of distribution
        ##

        ## -- domain 
        lb <- domain[1]
        ub <- domain[2]
        
        ## -- text for calling function
        cat(distr,"(",paste(...,sep=","),
            ") domain=(",signif(lb),",",signif(ub),"): ",sep="")

        ## -- function for generating a random sample
        if (is.null(rfunc)) {
                rfuncname <- paste("ur",distr,sep="")
                if (!exists(rfuncname))
                        stop("undefined function '",rfuncname,"'")
                rfunc <- match.fun(rfuncname, descend=FALSE)
        }
        ## -- function for computing probability vector
        if (is.null(dfunc) && missing(pv)) {
                dfuncname <- paste("d",distr,sep="")
                if (!exists(dfuncname))
                        stop("undefined function '",dfuncname,"'")
                dfunc <- match.fun(dfuncname, descend=FALSE)
        }
        
        ## -- run test and repeat once when failed the first time
        for (i in 1:2) {

                ## -- create probability vector
                if (missing(pv))
                        probs <- dfunc(lb:ub,...)
                else
                        probs <- pv

                ## -- random sample
                x <- rfunc(samplesize,lb=lb,ub=ub,...)

                ## -- test domain
                too.small <- length(x[x<lb])
                too.large <- length(x[x>ub])
                if (too.small > 1 || too.large > 1) {
                        too.small <- 100*too.small/samplesize
                        too.large <- 100*too.large/samplesize
                        stop("X out of domain (",too.small,"%|",too.large,"%)", call.=FALSE)
                }

                ## -- random distribution
                hits <- hist(x,breaks=(lb:(ub+1)-0.5),plot=FALSE)$counts

                ## -- collapse classes with too few entries
                expect <- samplesize * probs
                if (any (expect>5) ) {
                        accept <- which(expect > 5)
                        probs <- probs[accept]
                        hits <- hits[accept]
                        probs[1] <- probs[1] + sum(probs[-accept])
                        hits[1] <- hits[1] + sum(hits[-accept])
                }

                ## -- run unur.chiq.test --
                pval <- chisq.test(hits,p=probs,rescale.p=TRUE)$p.value

                ## -- check p-value
                if (pval > alpha) { # test passed
                        message("chisq test PASSed with p-value=",signif(pval))
                        break
                }

                ## -- test failed
                if (i>1) { # second try failed
                        stop("chisq test FAILED!  p-value=",signif(pval), call.=FALSE)
                }
                else {
                        warning("first chisq test FAILed with p-value=",signif(pval),
                                call.=FALSE,immediate.=TRUE)
                        assign("n.warns", get("n.warns",envir=unur.envir)+1, env=unur.envir)

                }
        }

        ## -- update list of p-values
        assign("pvals", append(get("pvals",envir=unur.envir), pval), env=unur.envir)

} ## --- end of unur.test.discr() ---


#############################################################################
##                                                                          #
##  Continuous multivariate Distributions                                   #
##                                                                          #
#############################################################################

## --- CMV: Function for running chi^2 goodness-of-fit test -----------------

unur.test.cmv <- function (distr, rfunc=NULL, pfunc=NULL, ...) {
        ##  Run a chi^2 test and evaluate p-value.
        ##  Repeat test once if it fails
        ##  (we do not check the validity of the algorithm
        ##   but the correct implementation.)
        ##  
        ##  Currently only the first component is tested!!!!
        ##
        ##  distr  ... name of distribution (as used for [rpqd]... calls
        ##  rfunc  ... random number generation function
        ##  pfunc  ... probability function (CDF) for marginal distributions
        ##  '...'  ... list of parameters of distribution
        ##

        ## -- text for calling function
        cat(distr,"(",paste(...,sep=","),")",
            sep="")

        ## -- function for generating a random sample
        if (is.null(rfunc))
                stop("sampling function required")
        ## -- function for computing CDF of marginal distributions
        if (is.null(pfunc))
                stop("probability function required")
        
        ## -- run test and repeat once when failed the first time
        for (i in 1:2) {

                ## -- random sample
                x <- rfunc(samplesize,...)
                ## -- transform into uniform distribution
                u <- pfunc(x[,1],...)
                ## -- make histogram of with equalsized bins (classified data)
                nbins <- as.integer(sqrt(samplesize))
                breaks <- (0:nbins)/nbins
                h <- hist(u,plot=F,breaks=breaks)$count
                ## -- run unur.chiq.test --
                pval <- chisq.test(h)$p.value
                ## -- check p-value
                if (pval > alpha) { # test passed
                        message("chisq test PASSed with p-value=",signif(pval))
                        break
                }

                ## -- test failed
                if (i>1) { # second try failed
                        stop("chisq test FAILED!  p-value=",signif(pval), call.=FALSE)
                }
                else {
                        warning("first chisq test FAILed with p-value=",signif(pval),
                                call.=FALSE,immediate.=TRUE)
                        assign("n.warns", get("n.warns",envir=unur.envir)+1, env=unur.envir)

                }
        }

        ## -- update list of p-values
        assign("pvals", append(get("pvals",envir=unur.envir), pval), env=unur.envir)

} ## --- end of unur.test.cmv() ---


#############################################################################
##                                                                          #
##  Auxiliary routines                                                      #
##                                                                          #
#############################################################################

## -- Print statistics ------------------------------------------------------

unur.test.statistic <- function () {
        pvals <- get("pvals", envir=unur.envir)
        n.warns <- get("n.warns", envir=unur.envir)
        cat("\nTests for discrete distributions\n\n",
            "\tnumber of tests = ",length(pvals),
            "(number of warnings = ",n.warns,")\n\n")
        summary(pvals)

        ## call garbage collector.
        ## this improves valgrind results on memory leaks
        silent <- gc()
}

## -- End -------------------------------------------------------------------
