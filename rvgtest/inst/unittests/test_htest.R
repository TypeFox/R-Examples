#############################################################################
##                                                                         ## 
## Test GoF tests with RVG frequency tables                                ## 
##                                                                         ## 
#############################################################################

## Test parameters ----------------------------------------------------------

## sample size
n <- 1e+4

## --------------------------------------------------------------------------
##
## Statistical tests for class 'rvgt.ftable'
##
## --------------------------------------------------------------------------

context("[htable] RVG frequency tables - GoF tests")

## rvgt.chisq ---------------------------------------------------------------

test_that("[ht-051] running rvgt.chisq: cont, pdist", {
        rep <- 5
        ft <- rvgt.ftable(n=n,rep=rep, rdist=rnorm,pdist=pnorm)
        ht <- rvgt.chisq(ft)
        pval <- ht$pval[rep]
        expect_true(pval > 1e-5, info="rvgt.chisq() broken")
})

## rvgt.Mtest ---------------------------------------------------------------

test_that("[ht-052] running rvgt.Mtest: cont, pdist", {
        rep <- 5
        ft <- rvgt.ftable(n=n,rep=rep, rdist=rnorm,pdist=pnorm)
        ht <- rvgt.Mtest(ft)
        pval <- ht$pval[rep]
        expect_true(pval > 1e-5, info="rvgt.Mtest() broken")
})

## plot.rvgt.htest ----------------------------------------------------------

test_that("[ht-061] plotting rvgt.htest", {
        ## we just run the code 
        ft <- rvgt.ftable(n=n,rep=5, rdist=rnorm,pdist=pnorm)
        ht1 <- rvgt.chisq(ft)
        retval <- plot(ht1)
        expect_identical(retval,NULL)
        
        ht2 <- rvgt.Mtest(ft)
        retval <- plot(ht2)
        expect_identical(retval,NULL)
        
        retval <- plot.rvgt.htest(list(ht1,ht2))
        expect_identical(retval,NULL)
        
        ft <- rvgt.ftable(n=n,rep=1, rdist=rnorm,pdist=pnorm)
        ht3 <- rvgt.chisq(ft)
        retval <- plot(ht3)
        expect_identical(retval,NULL)
        
        retval <- plot.rvgt.htest(list(ht1,ht2,ht3))
        expect_identical(retval,NULL)
})


## --------------------------------------------------------------------------
##
## Check invalid arguments
##
## --------------------------------------------------------------------------

context("[htable] RVG frequency tables - Invalid arguments for GoF tests")

## rvgt.chisq ---------------------------------------------------------------

test_that("[ht-i21] calling rvgt.chisq with invalid arguments", {
        msg <- "Argument 'ftable' missing or invalid."
        expect_error(rvgt.chisq(),          msg)
        expect_error(rvgt.chisq("ftable"),  msg)
        expect_error(rvgt.chisq(list(a=1)), msg)
})

## rvgt.Mtest ---------------------------------------------------------------

test_that("[ht-i22] calling rvgt.Mtest with invalid arguments", {
        msg <- "Argument 'ftable' missing or invalid."
        expect_error(rvgt.Mtest(),          msg)
        expect_error(rvgt.Mtest("ftable"),  msg)
        expect_error(rvgt.Mtest(list(a=1)), msg)
})

## plot.rvgt.htest ----------------------------------------------------------

test_that("[ht-i41] calling plot.rvgt.htest with invalid arguments: alpha", {
        ## significance level 'alpha'
        ft <- rvgt.ftable(n=n, rdist=rnorm, qdist=qnorm)
        ht <- rvgt.chisq(ft)

        msg <- "Invalid argument 'alpha'."
        expect_error(plot(ht,alpha="a"), msg)
        expect_error(plot(ht,alpha=0),   msg)
        expect_error(plot(ht,alpha=1),   msg)
})

test_that("[ht-i42] calling plot.rvgt.htest with invalid arguments: x", {
        ## htest object
        ft <- rvgt.ftable(n=n, rdist=rnorm, qdist=qnorm)
        ht <- rvgt.chisq(ft)

        msg <- "Argument 'x' missing or invalid."
        expect_error(plot.rvgt.htest(               alpha=0.05), msg)
        expect_error(plot.rvgt.htest(x=1:10,        alpha=0.05), msg)

        msg <- "Invalid argument 'x'."
        expect_error(plot.rvgt.htest(x=list(a=1:10),alpha=0.05), msg)
})


## -- Clear -----------------------------------------------------------------

rm(n)

## -- End -------------------------------------------------------------------
