#############################################################################
##                                                                         ## 
## Test functions for estimating error in numerical inversion methods      ##
##                                                                         ## 
#############################################################################

## Test parameters ----------------------------------------------------------

## Auxiliary routines -------------------------------------------------------

check.uerror <- function (error, limit=1e-10) {
        error_max <- max(error$max)
        expect_true(error_max < limit, info="uerror() broken", label="error < limit")
}

check.xerror <- function (error, limit=1e-10) {
        error_max <- max(error$max)
        expect_true(error_max < limit, info="xerror() broken", label="error < limit")
}


## --------------------------------------------------------------------------
##
## uerror()
##
## --------------------------------------------------------------------------

context("[ierror] Approximation error - uerror")

## uerror -------------------------------------------------------------------

test_that("[ie-001] calling uerror: default", {
        ue <- uerror(n=1e4, aqdist=qnorm, pdist=pnorm)
        check.uerror(ue)
})

test_that("[ie-002] calling uerror: res", {
        ue <- uerror(n=1e4, res=100, aqdist=qnorm, pdist=pnorm)
        check.uerror(ue)
})

test_that("[ie-003] calling uerror: tails", {
        ue <- uerror(n=1e3, res=100, aqdist=qnorm, pdist=pnorm, tails=TRUE)
        check.uerror(ue)
})

test_that("[ie-004] calling uerror: udomain", {
        ue <- uerror(n=1e3, res=100, aqdist=qnorm, pdist=pnorm, udomain=c(0,0.1))
        check.uerror(ue)
})

test_that("[ie-005] calling uerror: dparams", {
        ue <- uerror(n=1e3, res=100, aqdist=function(u){qgamma(u,shape=2)},
                     pdist=pgamma, shape=2, udomain=c(0,0.1))
        check.uerror(ue)
})

test_that("[ie-006] calling uerror: dparams, plot", {
        ue <- uerror(n=1e3, res=100, aqdist=function(u){qgamma(u,shape=2)},
                     pdist=pgamma, shape=2, udomain=c(0,0.1), plot=TRUE)
        check.uerror(ue)
})

test_that("[ie-007] calling uerror: trunc", {
        ## An inverse CDF for a truncated normal distribution
        aqtn <- function(x) { qnorm(x * (pnorm(2.5) - pnorm(1.5)) + pnorm(1.5)) }

        ue <- uerror(n=1e5, res=100, aqdist=aqtn, pdist=pnorm, trunc=c(1.5,2.5))
        check.uerror(ue)
})

## --------------------------------------------------------------------------
##
## xerror()
##
## --------------------------------------------------------------------------

context("[ierror] Approximation error - xerror")

## xerror -------------------------------------------------------------------

test_that("[ie-011] calling xerror: default", {
        aq <- function(u) { qnorm(u) + u*runif(length(u),max=9.e-11) }
        xe <- xerror(n=1e5, aqdist=aq, qdist=qnorm)
        check.xerror(xe)
})

test_that("[ie-012] calling xerror: res", {
        aq <- function(u) { qnorm(u) + u*runif(length(u),max=9.e-11) }
        xe <- xerror(n=1e4, res=100, aqdist=aq, qdist=qnorm)
        check.xerror(xe)
})

test_that("[ie-013] calling xerror: tails", {
        aq <- function(u) { qnorm(u) + u*runif(length(u),max=9.e-11) }
        xe <- xerror(n=1e3, res=100, aqdist=aq, qdist=qnorm, tails=TRUE)
        check.xerror(xe)
})

test_that("[ie-014] calling xerror: udomain", {
        aq <- function(u) { qnorm(u) + u*runif(length(u),max=9.e-11) }
        xe <- xerror(n=1e4, res=100, aqdist=aq, qdist=qnorm, udomain=c(0,0.01))
        check.xerror(xe)
})

test_that("[ie-015] calling xerror: plot", {
        aq <- function(u) { qnorm(u) + u*runif(length(u),max=9.e-11) }
        xe <- xerror(n=1e4, res=100, aqdist=aq, qdist=qnorm, plot=TRUE)
        check.xerror(xe)
})

test_that("[ie-016] calling xerror: abs", {
        aq <- function(u) { qnorm(u) + u*runif(length(u),max=9.e-11) }
        xe <- xerror(n=1e4, res=100, kind="abs", aqdist=aq, qdist=qnorm, plot=TRUE)
        check.xerror(xe)
})

test_that("[ie-017] calling xerror: rel", {
        aq <- function(u) { qnorm(u) + u*runif(length(u),max=9.e-11) }
        xe <- xerror(n=1e4, res=100, kind="rel", aqdist=aq, qdist=qnorm, plot=TRUE)
        check.xerror(xe, limit=1e-5)
})

test_that("[ie-018] calling xerror: abs, trunc", {
        aqtn <- function(x) { qnorm(x * (pnorm(2.5) - pnorm(1.5)) + pnorm(1.5)) }
        xe <- xerror(n=1e5, res=100, aqdist=aqtn, qdist=qnorm, trunc=c(1.5,2.5), kind="abs")
        check.xerror(xe)
})

test_that("[ie-019] calling xerror: rel, trunc", {
  aqtn <- function(x) { qnorm(x * (pnorm(2.5) - pnorm(1.5)) + pnorm(1.5)) }
  xe <- xerror(n=1e5, res=100, aqdist=aqtn, qdist=qnorm, trunc=c(1.5,2.5), kind="rel")
  check.xerror(xe)
})


## --------------------------------------------------------------------------
##
## plot.rvgt.ierror()
##
## --------------------------------------------------------------------------

context("[ierror] Approximation error - plot")

## plot.rvgt.ierror ---------------------------------------------------------

test_that("[ie-021] calling plot.rvgt.ierror: uerror", {
        ## we just run the code 
        ue1 <- uerror(n=1e3, res=100, aqdist=function(u){qgamma(u,shape=2)},
                      pdist=pgamma, shape=2)
        retval <- plot(ue1)
        expect_identical(retval,NULL)

        retval <- plot(ue1,maxonly=TRUE)
        expect_identical(retval,NULL)
  
        ue2 <- uerror(n=1e3, res=100, aqdist=function(u){qgamma(u,shape=1)},
                      pdist=pgamma, shape=1, udomain=c(0.1,0.8))

        retval <- plot(ue2)
        expect_identical(retval,NULL)
  
        retval <- plot(ue2,tol=1e-15)
        expect_identical(retval,NULL)
  
        retval <- plot(ue2,tol=1e-16)
        expect_identical(retval,NULL)
  
        retval <- plot(ue2,tol=5e-18)
        expect_identical(retval,NULL)
  
        retval <- plot.rvgt.ierror(list(ue1,ue2),tol=4.e-16)
        expect_identical(retval,NULL)
})

test_that("[ie-022] calling plot.rvgt.ierror: xerror, abs", {
        ## we just run the code 
        aq1 <- function(u) { qnorm(u) + u*runif(length(u),max=1.e-10) }
        xe1 <- xerror(n=1e3, res=100, aqdist=aq1, qdist=qnorm)

        retval <- plot(xe1)
        expect_identical(retval,NULL)

        retval <- plot(xe1,maxonly=TRUE)
        expect_identical(retval,NULL)
        
        aq2 <- function(u) { qnorm(u) + (1-u*u)*runif(length(u),max=1.e-10) }
        xe2 <- xerror(n=1e3, res=100, aqdist=aq2, qdist=qnorm)

        retval <- plot(xe2)
        expect_identical(retval,NULL)

        retval <- plot.rvgt.ierror(list(xe1,xe2),tol=5.e-11)
        expect_identical(retval,NULL)
})

test_that("[ie-023] calling plot.rvgt.ierror: xerror, rel", {
        ## we just run the code 
        aq1 <- function(u) { qnorm(u) + u*runif(length(u),max=1.e-10) }
        xe1 <- xerror(n=1e3, res=100, aqdist=aq1, qdist=qnorm, kind="rel")
        retval <- plot(xe1)
        expect_identical(retval,NULL)
        
        aq2 <- function(u) { qnorm(u) + (1-u*u)*runif(length(u),max=1.e-10) }
        xe2 <- xerror(n=1e3, res=100, aqdist=aq2, qdist=qnorm, kind="rel")
        retval <- plot(xe2)
        expect_identical(retval,NULL)
        
        retval <- plot.rvgt.ierror(list(xe1,xe2),tol=5.e-11)
        expect_identical(retval,NULL)
})


## --------------------------------------------------------------------------
##
## Check invalid arguments
##
## --------------------------------------------------------------------------

context("[ierror] Approximation error - Invalid arguments")

## uerror -------------------------------------------------------------------

test_that("[ie-i11] calling uerror with invalid arguments: n", {
        ## sample size 'n'
        msg <- "Argument 'n' missing or invalid."
        expect_error(uerror(       aqdist=qnorm, pdist=pnorm), msg)
        expect_error(uerror(n="a", aqdist=qnorm, pdist=pnorm), msg)
        expect_error(uerror(n=0  , aqdist=qnorm, pdist=pnorm), msg)
        expect_error(uerror(n=1.2, aqdist=qnorm, pdist=pnorm), msg)
})

test_that("[ie-i12] calling uerror with invalid arguments: res", {
        ## resolution 'res'
        msg <- "Invalid argument 'res'."
        expect_error(uerror(n=1e4, res="a", aqdist=qnorm, pdist=pnorm), msg)
        expect_error(uerror(n=1e4, res=0,   aqdist=qnorm, pdist=pnorm), msg)
        expect_error(uerror(n=1e4, res=1.2, aqdist=qnorm, pdist=pnorm), msg)

        msg <- "Invalid arguments 'n' < 'res'."
        expect_error(uerror(n=100, res=201, aqdist=qnorm, pdist=pnorm), msg)

})

test_that("[ie-i13] calling uerror with invalid arguments: aqdist", {
        ## approximate inverse distribution function (quantile function)
        msg <- "Argument 'aqdist' missing."
        expect_error(uerror(n=1e3, res=100                 ), msg)

        msg <- "Argument 'aqdist' invalid."
        expect_error(uerror(n=1e3, res=100, aqdist="aqdist"), msg)
})

test_that("[ie-i14] calling uerror with invalid arguments: pdist", {
        ## distribution function
        msg <- "Argument 'pdist' missing or invalid."
        expect_error(uerror(n=1e3, res=100, aqdist=qnorm               ), msg)
        expect_error(uerror(n=1e3, res=100, aqdist=qnorm, pdist="pnorm"), msg)
})

test_that("[ie-i15] calling uerror with invalid arguments: udomain", {
        ## domain
        msg <- "Invalid argument 'udomain'."
        expect_error(uerror(n=1e3, res=100, aqdist=qnorm, pdist=pnorm, udomain=c(0.5,0.5)), msg)
        expect_error(uerror(n=1e3, res=100, aqdist=qnorm, pdist=pnorm, udomain=c(-0.5,0.)), msg)
})


## xerror -------------------------------------------------------------------

test_that("[ie-i21] calling xerror with invalid arguments: n", {
        ## sample size 'n'
        msg <- "Invalid argument 'n'."
        expect_error(xerror(       aqdist=qnorm, qdist=qnorm), msg)
        expect_error(xerror(n="a", aqdist=qnorm, qdist=qnorm), msg)
        expect_error(xerror(n=0  , aqdist=qnorm, qdist=qnorm), msg)
        expect_error(xerror(n=1.2, aqdist=qnorm, qdist=qnorm), msg)
})

test_that("[ie-i22] calling xerror with invalid arguments: res", {
        ## resolution 'res'
        msg <- "Invalid argument 'res'."
        expect_error(xerror(n=1e4, res="a", aqdist=qnorm, qdist=qnorm), msg)
        expect_error(xerror(n=1e4, res=0,   aqdist=qnorm, qdist=qnorm), msg)
        expect_error(xerror(n=1e4, res=1.2, aqdist=qnorm, qdist=qnorm), msg)
        msg <- "Invalid arguments 'n' < 'res'."
        expect_error(xerror(n=100, res=201, aqdist=qnorm, qdist=qnorm), msg)
})

test_that("[ie-i23] calling xerror with invalid arguments: aqdist", {
        ## approximate inverse distribution function (quantile function)
        msg <- "Argument 'aqdist' missing."
        expect_error(xerror(n=1e3, res=100                 ), msg)
        msg <- "Argument 'aqdist' invalid."
        expect_error(xerror(n=1e3, res=100, aqdist="aqdist"), msg)
})

test_that("[ie-i24] calling xerror with invalid arguments: qdist", {
        ## (exact) quatile function of distribution
        msg <- "Argument 'qdist' missing or invalid."
        expect_error(xerror(n=1e3, res=100, aqdist=qnorm               ), msg)
        expect_error(xerror(n=1e3, res=100, aqdist=qnorm, qdist="qnorm"), msg)
})

test_that("[ie-i25] calling xerror with invalid arguments: udomain", {
        ## domain
        msg <- "Invalid argument 'udomain'."
        expect_error(xerror(n=1e3, res=100, aqdist=qnorm, qdist=qnorm, udomain=c(0.5,0.5)), msg)
        expect_error(xerror(n=1e3, res=100, aqdist=qnorm, qdist=qnorm, udomain=c(-0.5,0.)), msg)
})

test_that("[ie-i26] calling xerror with invalid arguments: kind", {
        ## kind
        msg <- "'arg' should be one of \"abs\", \"rel\""
        expect_error(xerror(n=1e4, aqdist=qnorm, qdist=qnorm, kind="foo"), msg)
        msg <- "'arg' must be NULL or a character vector"
        expect_error(xerror(n=1e4, aqdist=qnorm, qdist=qnorm, kind=10000), msg)
})

## plot.rvgt.ierror ---------------------------------------------------------

test_that("[ie-i31] calling plot.rvgt.ierror: x", {
        msg <- "'x' is missing."
        expect_error(plot.rvgt.ierror(              ), msg)

        msg <- "Invalid argument 'x'."
        expect_error(plot.rvgt.ierror(x=1:10        ), msg)
        expect_error(plot.rvgt.ierror(x=list(a=1:10)), msg)
})

test_that("[ie-i32] calling plot.rvgt.ierror: tol", {
        msg <- "Invalid argument 'tol'."
        ue <- uerror(n=1e4, res=10, aqdist=qnorm, pdist=pnorm)
        expect_error(plot.rvgt.ierror(ue, tol=0   ), msg)
        expect_error(plot.rvgt.ierror(ue, tol=-0.1), msg)
})

## -- End -------------------------------------------------------------------
