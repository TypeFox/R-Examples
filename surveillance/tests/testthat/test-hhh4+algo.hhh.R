context("Comparison of hhh4() and algo.hhh() for 'influMen' example")

## influenza/meningococcal data (also illustrated in vignette("hhh4")
data("influMen")

## fit with old algo.hhh()
hhhfit <- algo.hhh(influMen,
                   list(lambda=c(1,1), neighbours=c(NA,0),
                        linear=FALSE, nseason=c(3,1), negbin="multiple"),
                   verbose=FALSE)
test_that("algo.hhh() converges for 'influMen' example",
          expect_true(hhhfit$convergence))

## fit with new hhh4()
hhh4fit <- hhh4(disProg2sts(influMen),
                list(ar=list(f=~0+fe(1, which=c(TRUE, TRUE)), lag=1),
                     ne=list(f=~0+fe(1, which=c(FALSE,TRUE)), lag=0,
                             weights=matrix(c(0,0,1,0), 2, 2)), # influenza->IMD
                     end=list(f=addSeason2formula(
                              ~0+fe(1,which=c(TRUE,TRUE)), S=c(3,1))),
                     family="NegBinM"))
test_that("hhh4() converges for 'influMen' example",
          expect_true(hhh4fit$convergence))

## compare fits
test_that("results from algo.hhh() and hhh4() agree for 'influMen' example", {
    expect_that(hhhfit$nObs, equals(nobs(hhh4fit)))
    expect_that(hhhfit$loglikelihood, equals(hhh4fit$loglikelihood))
    ## fitted values
    expect_that(fitted(hhhfit),
                equals(fitted(hhh4fit),
                       tolerance=0.0005, check.attributes=FALSE))
    ## coefficient estimates
    hhh4coefs <- coef(hhh4fit, idx2Exp=1:3, reparamPsi=TRUE)
    orderhhh42old <- c(4,5,1:3,
                       grep("(sin|cos).*\\.influenza", names(hhh4coefs)),
                       grep("(sin|cos).*\\.meningo", names(hhh4coefs)),
                       grep("overdisp", names(hhh4coefs)))
    expect_that(coef(hhhfit, reparamPsi=TRUE),
                equals(hhh4coefs[orderhhh42old],
                       tolerance=0.0005, check.attributes=FALSE))
    ## variance-covariance matrix of parameter estimates
    hhh4cov <- vcov(hhh4fit, idx2Exp=c(1:3,14:15), reparamPsi=FALSE)
    expect_that(hhhfit$cov,
                equals(hhh4cov[orderhhh42old,orderhhh42old],
                       tolerance=0.002, check.attributes=FALSE))
})
