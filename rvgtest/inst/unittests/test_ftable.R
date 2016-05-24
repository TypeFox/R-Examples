#############################################################################
##                                                                         ## 
## Test RVG frequency table                                                ## 
##                                                                         ## 
#############################################################################

## Test parameters ----------------------------------------------------------

## sample size
n <- 1e+4

## Auxiliary routines -------------------------------------------------------

check.ftable <- function (ftable, rep=1) {
        pval <- rvgt.chisq(ftable)$pval[rep]
        expect_true(pval > 1e-5, info="rvgt.ftable() broken")
}
 
## --------------------------------------------------------------------------
##
## rvgt.ftable()
##
## --------------------------------------------------------------------------

context("[ftable] RVG frequency tables - Create")

## rvgt.ftable --------------------------------------------------------------

test_that("[ft-001] calling rvgt.ftable: default", {
        ft <- rvgt.ftable(n=n, rdist=rnorm,qdist=qnorm)
        check.ftable(ft)
})

test_that("[ft-002] calling rvgt.ftable: dparams", {
        ft <- rvgt.ftable(n=n, rdist=rnorm,qdist=qnorm, mean=1,sd=2)
        check.ftable(ft)
})

test_that("[ft-003] calling rvgt.ftable: dparams, breaks", {
        ft <- rvgt.ftable(n=n, rdist=rnorm,qdist=qnorm, breaks=51, mean=1,sd=2)
        check.ftable(ft)
        ft <- rvgt.ftable(n=n, rdist=rnorm,qdist=qnorm, breaks=(0:100)/100)
        check.ftable(ft)
})

test_that("[ft-004] calling rvgt.ftable: dparams, breaks, rep", {
        ft <- rvgt.ftable(n=n, rep=5, rdist=rnorm,qdist=qnorm, breaks=51, mean=1,sd=2)
        check.ftable(ft, rep=5)
})

test_that("[ft-005] calling rvgt.ftable: exactu", {
        ft <- rvgt.ftable(n=n,rep=1, rdist=rnorm,pdist=pnorm, exactu=TRUE)
        check.ftable(ft)
        ft <- rvgt.ftable(n=n,rep=1, rdist=rnorm,pdist=pnorm, exactu=FALSE)
        check.ftable(ft)
})

test_that("[ft-006] calling rvgt.ftable: plot", {
        ft <- rvgt.ftable(n=n,rep=5, rdist=rnorm,qdist=qnorm, breaks=51, plot=TRUE)
        check.ftable(ft, rep=5)
})

## univariate continuous distribution ---------------------------------------
 
test_that("[ft-011] calling rvgt.ftable: cont, pdist", {
        ft <- rvgt.ftable(n=n, rdist=rnorm,pdist=pnorm)
        check.ftable(ft)
})

test_that("[ft-012] calling rvgt.ftable: cont, pdist, exactu", {
        ft <- rvgt.ftable(n=n, rdist=rnorm,pdist=pnorm, exactu=TRUE)
        check.ftable(ft)
})

test_that("[ft-013] calling rvgt.ftable: cont, qdist", {
        ft <- rvgt.ftable(n=n, rdist=rnorm,qdist=qnorm)
        check.ftable(ft)
})

test_that("[ft-014] calling rvgt.ftable: cont, pdist, qdist", {
   ft <- rvgt.ftable(n=n, rdist=rnorm,qdist=qnorm,pdist=pnorm)
   check.ftable(ft)
})

## univariate discrete distribution -----------------------------------------

test_that("[ft-021] calling rvgt.ftable: discr, pdist, dparams", {
        ft <- rvgt.ftable(n=n, rdist=rgeom,pdist=pgeom, prob=0.123)
        check.ftable(ft)
})

## 'qdist' without 'pdist' does not work for discrete distributions.

## truncated domain ---------------------------------------------------------

test_that("[ft-031] calling rvgt.ftable: cont, trunc, pdist, qdist", {
        rdist <- function(n) {
                x <- numeric(n)
                for (i in 1:n) {
                        while(TRUE) { x[i] <- rnorm(1); if (x[i]>0 && x[i]<1) break } }
                return(x)
        }
 
        ft <- rvgt.ftable(n=1e3,rep=5, rdist=rdist, pdist=pnorm, qdist=qnorm, trunc=c(0,1))
        check.ftable(ft, rep=5)
})

test_that("[ft-032] calling rvgt.ftable: cont, trunc, qdist", {
        rdist <- function(n) {
                x <- numeric(n)
                for (i in 1:n) {
                        while(TRUE) { x[i] <- rnorm(1); if (x[i]>0 && x[i]<1) break } }
                return(x)
        }
 
        ft <- rvgt.ftable(n=1e3,rep=5, rdist=rdist, qdist=qnorm, trunc=c(0,1))
        check.ftable(ft, rep=5)
})

test_that("[ft-033] calling rvgt.ftable: cont, trunc, pdist", {
        rdist <- function(n) {
                x <- numeric(n)
                for (i in 1:n) {
                        while(TRUE) { x[i] <- rnorm(1); if (x[i]>0 && x[i]<1) break } }
                return(x)
        }
 
        ft <- rvgt.ftable(n=2e3,rep=5, rdist=rdist, pdist=pnorm, trunc=c(0,1))
        check.ftable(ft, rep=5)
})

## There currently no tests for truncated discrete distributions.

## plot.rvgt.ftable ---------------------------------------------------------

test_that("[ft-041] plotting rvgt.ftable", {
        ## we just run the code 
        ft <- rvgt.ftable(n=n,rep=5, rdist=rnorm,pdist=pnorm, exactu=TRUE)
        retval <- plot(ft)
        expect_identical(retval,NULL)
        
        retval <- plot(ft,rows=c(2,3),alpha=0.005)
        expect_identical(retval,NULL)

        ft <- rvgt.ftable(n=n,rep=5, rdist=rnorm,pdist=pnorm, exactu=FALSE)
        retval <- plot(ft)
        expect_identical(retval,NULL)
})


## --------------------------------------------------------------------------
##
## Check invalid arguments
##
## --------------------------------------------------------------------------

context("[ftable] RVG frequency tables - Invalid arguments")

## rvgt.ftable --------------------------------------------------------------
 
test_that("[ft-i11] calling rvgt.ftable with invalid arguments: n", {
        ## sample size 'n'
        msg <- "Argument 'n' missing or invalid."
        expect_error(rvgt.ftable(       rdist=rnorm, qdist=qnorm), msg)
        expect_error(rvgt.ftable(n="a", rdist=rnorm, qdist=qnorm), msg)
        expect_error(rvgt.ftable(n=0,   rdist=rnorm, qdist=qnorm), msg)
        expect_error(rvgt.ftable(n=1.2, rdist=rnorm, qdist=qnorm), msg)
})

test_that("[ft-i12] calling rvgt.ftable with invalid arguments: res", {
        ## resolution 'res'
        msg <- "Invalid argument 'rep'."
        expect_error(rvgt.ftable(n=100, rep="a", rdist=rnorm, qdist=qnorm), msg)
        expect_error(rvgt.ftable(n=100, rep=0,   rdist=rnorm, qdist=qnorm), msg)
        expect_error(rvgt.ftable(n=100, rep=1.2, rdist=rnorm, qdist=qnorm), msg)
})

test_that("[ft-i13] calling rvgt.ftable with invalid arguments: rdist", {
        ## random variate generator
        msg <- "Argument 'rdist' missing."
        expect_error(rvgt.ftable(n=100, qdist=rnorm               ), msg)

        msg <- "Argument 'rdist' invalid."
        expect_error(rvgt.ftable(n=100, rdist="rnorm", qdist=qnorm), msg)
})

test_that("[ft-i14] calling rvgt.ftable with invalid arguments: pdist, qdist", {
        ## quantile and distribution function: continuous
        msg <- "Argument 'qdist' or 'pdist' required."
        expect_error(rvgt.ftable(n=100, rdist=rnorm               ), msg)

        msg <- "Argument 'pdist' invalid."
        expect_error(rvgt.ftable(n=100, rdist=rnorm, pdist="pnorm"), msg)

        msg <- "Argument 'qdist' invalid."
        expect_error(rvgt.ftable(n=100, rdist=rnorm, qdist="qnorm"), msg)

        ## quantile and distribution function: discrete
        msg <- "Argument 'pdist' required for discrete distribution."
        expect_error(rvgt.ftable(n=100, rdist=rgeom,qdist=qgeom, prob=0.123), msg)

        msg <- "Argument 'qdist' ignored for discrete distributions."
        expect_warning_suppress(
          rvgt.ftable(n=n, rdist=rgeom,pdist=pgeom,qdist=qgeom, prob=0.123), msg)
})

test_that("[ft-i15] calling rvgt.ftable with invalid arguments: breaks", {
        ## break points
        msg <- "Invalid argument 'breaks'."
        expect_error(rvgt.ftable(n=100, rdist=rnorm, qdist=qnorm, breaks=c(0,0.1,0.2,"0.3")), msg)
        expect_error(rvgt.ftable(n=100, rdist=rnorm, qdist=qnorm, breaks=numeric()),          msg)

        msg <- "Number of break points too small"
        expect_error(rvgt.ftable(n=100, rdist=rnorm, qdist=qnorm, breaks=2),                  msg)

        msg <- "Number of break points too small"
        expect_error(rvgt.ftable(n=100, rdist=rnorm, qdist=qnorm, breaks=c(0,1)),             msg)

        msg <- "break points out of range"
        expect_error(rvgt.ftable(n=100, rdist=rnorm, qdist=qnorm, breaks=c(0,-2,1)),          msg)

        msg <- "break points invalid: length of histogram cells"
        expect_error(rvgt.ftable(n=100, rdist=rnorm, qdist=qnorm, breaks=c(0,0.5,0.5,1)),     msg)

        ## use exact location of break points
        msg <- "Argument 'exactu' must be boolean."
        expect_error(rvgt.ftable(n=100, rdist=rnorm, qdist=qnorm, exactu=0), msg)
})

test_that("[ft-i16] calling rvgt.ftable with invalid arguments: exactu", {
        ## exactu
        msg <- "Argument 'exactu' must be boolean."
        expect_error(rvgt.ftable(n=100, rdist=rnorm, qdist=qnorm, exactu="on"), msg)

        msg <- "Argument 'exactu' ignored for discrete distributions."
        expect_warning_suppress(
          rvgt.ftable(n=n, rdist=rgeom,pdist=pgeom,exactu=TRUE, prob=0.123), msg)
})

## plot.rvgt.ftable ---------------------------------------------------------

test_that("[ft-i31] calling plot.rvgt.ftable with invalid arguments: rows", {
        ## number of rows
        ft <- rvgt.ftable(n=n, rdist=rnorm, qdist=qnorm)
        
        msg <- "Invalid argument 'rows'."
        expect_error(plot(ft,rows="a"), msg)
        expect_error(plot(ft,rows=0),   msg)
        expect_error(plot(ft,rows=2),   msg)
})

test_that("[ft-i32] calling plot.rvgt.ftable with invalid arguments: alpha", {
        ## significance level 'alpha'
        ft <- rvgt.ftable(n=n, rdist=rnorm, qdist=qnorm)

        msg <- "Invalid argument 'alpha'."
        expect_error(plot(ft,alpha="a"), msg)
        expect_error(plot(ft,alpha=0),   msg)
        expect_error(plot(ft,alpha=1),   msg)
})

## -- Clear -----------------------------------------------------------------

rm(n)

## -- End -------------------------------------------------------------------
