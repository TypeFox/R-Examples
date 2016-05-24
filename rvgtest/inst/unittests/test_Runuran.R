#############################################################################
##                                                                         ## 
## Test RVG interface for UNU.RAN objects                                  ## 
##                                                                         ## 
#############################################################################

## Load library -------------------------------------------------------------

if (require(Runuran)) {

## --------------------------------------------------------------------------
##
## Test functions
##
## --------------------------------------------------------------------------

context("[Runuran] Runuran API")

## rvgt.ftable --------------------------------------------------------------

check.ftable <- function (ftable, rep=1) {
        pval <- rvgt.chisq(ftable)$pval[rep]
        expect_true(pval > 1e-5, info="rvgt.ftable(): Runuran API does not work")
}

test_that("[Ru-001] calling rvgt.table with Runuran: cont, rdist", {
        unr <- new("unuran", "normal(1,2)")
        ft <- rvgt.ftable(n=1e5,rep=5, rdist=unr,qdist=qnorm, breaks=51, mean=1,sd=2)
        check.ftable(ft)
})

test_that("[Ru-002] calling rvgt.table with Runuran: discr, rdist", {
        unr <- dgtd.new(udbinom(size=20,prob=0.3))
        ft <- rvgt.ftable(n=1e5,rep=5, rdist=unr,pdist=pbinom, breaks=51, size=20, prob=0.3)
        check.ftable(ft)
})

test_that("[Ru-003] calling rvgt.table with Runuran: cont, rdist, trunc", {
        ## truncated domain
        unr <- new("unuran", "normal();domain=(0,1)")
        ft <- rvgt.ftable(n=1e5,rep=5, rdist=unr,qdist=qnorm, breaks=51, trunc=c(0,1))
        check.ftable(ft)
})

## uerror -------------------------------------------------------------------

check.uerror <- function (error, limit=1e-10) {
        error_max <- max(error$max)
        expect_true(error_max < limit, info="uerror(): Runuran API does not work", label="error < limit")
}

test_that("[Ru-011] calling uerror with Runuran: cont, aqdist", {
  unr <- pinvd.new(udnorm(mean=1, sd=2), uresolution=1.e-12)
  ue <- uerror(n=1e3, res=100, aqdist=unr, pdist=pnorm, mean=1,sd=2)
  check.uerror(ue)
})

## xerror -------------------------------------------------------------------

check.xerror <- function (error, limit=1e-10) {
        error_max <- max(error$max)
        expect_true(error_max < limit, info="xerror(): Runuran API does not work", label="error < limit")
}
          
test_that("[Ru-021] calling xerror with Runuran: cont, aqdist", {
        unr <- pinvd.new(udnorm(mean=1,sd=2), uresolution=1.e-12)
        xe <- xerror(n=1e3, res=100, aqdist=unr, qdist=qnorm, mean=1,sd=2)
        check.xerror(xe, limit=1e-8)
})


## --------------------------------------------------------------------------
##
## Check invalid arguments
##
## --------------------------------------------------------------------------

context("[Runuran] Runuran API - Invalid arguments")

## UNU.RAN object must be of class "unuran.cont" ----------------------------

test_that("[Ru-i11] calling rvgt.ftable (Runuran) with invalid arguments: multivariate", {
        ## rvgt.ftable() requires univariate distribution
        unr <- hitro.new(dim=2, pdf=function(x){exp(-sum(x^2))})
        msg <- "Argument 'rdist' is object of class 'unuran' of invalid distribution type."
        expect_error(rvgt.ftable(n=1e4,rep=1, rdist=unr,qdist=qbinom, size=20,prob=0.3), msg)
})

test_that("[Ru-i12] calling uerror (Runuran) with invalid arguments: discrete", {
        ## uerror() requires univariate continuous distribution
        unr <- dgtd.new(udbinom(size=20,prob=0.3))
        msg <- " Argument 'aqdist' is object of class 'unuran' of invalid distribution type."
        expect_error(uerror(n=1e3,res=100, aqdist=unr,pdist=pbinom, size=20,prob=0.3), msg)
})

test_that("[Ru-i12] calling xerror (Runuran) with invalid arguments: discrete", {
        ## xerror() requires univariate continuous distribution
        unr <- dgtd.new(udbinom(size=20,prob=0.3))
        msg <- " Argument 'aqdist' is object of class 'unuran' of invalid distribution type."
        expect_error(xerror(n=1e3,res=100, aqdist=unr,qdist=qbinom, size=20,prob=0.3), msg)
})

## UNU.RAN object with non-inversion method ---------------------------------

test_that("[Ru-i21] calling uerror (Runuran) with invalid arguments: TDR", {
        ## uerror() requires inversion method but TDR is rejection method
        unr <- tdrd.new(udnorm())
        msg <- "Argument 'aqdist' is invalid UNU.RAN object: inversion method required."
        expect_error(uerror(n=1e3, res=100, aqdist=unr, pdist=pnorm), msg)
})

test_that("[Ru-i22] calling xerror (Runuran) with invalid arguments: TDR", {
        ## xerror() requires inversion method but TDR is rejection method
        unr <- tdrd.new(udnorm())
        msg <- "Argument 'aqdist' is invalid UNU.RAN object: inversion method required."
        expect_error(xerror(n=1e3, res=100, aqdist=unr, qdist=qnorm), msg)
})


## --------------------------------------------------------------------------
##
## Package not available
##
## --------------------------------------------------------------------------

} else {
  message("\nCannot run test file 'test_Runuran' -- package 'Runuran' is not available\n")
} ## end of 'if (require(Runuran)) {'

## -- End -------------------------------------------------------------------
