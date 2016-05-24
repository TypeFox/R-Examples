context("Fixed effects hhh4() model fit and involved analytical derivatives")

data("measlesWeserEms")
measlesModel <- list(
    end = list(f = addSeason2formula(~1 + t, S=1, period=52),
               offset = population(measlesWeserEms)),
    ar = list(f = ~1),
    ne = list(f = ~1 + log(pop),
        weights = W_powerlaw(maxlag = 5, normalize = TRUE)),
    family = "NegBin1", data = list(pop = population(measlesWeserEms))
)

measlesFit <- hhh4(stsObj = measlesWeserEms, control = measlesModel)

test_that("estimates and standard errors are reproducible", {
    ## dput(coef(measlesFit, se = TRUE))
    orig <- structure(
        c(-0.499636482022272, 0.551345030080107, 0.96093157194767, 
          -0.153585641356373, 0.00333284018297979, 1.01500011496702,
          -0.588738943313705, 5.52782609236691, 1.81915612994789,
          0.121781347106564, 1.27401298230559, 0.453889365025671,
          0.281013375484401, 0.00459840327748742, 0.210642721317572, 
          0.191921649336323, 1.87984346848385, 0.265016986696184),
        .Dim = c(9L, 2L),
        .Dimnames = list(c("ar.1", "ne.1", "ne.log(pop)", "end.1", 
            "end.t", "end.sin(2 * pi * t/52)", "end.cos(2 * pi * t/52)", 
            "neweights.d", "overdisp"), c("Estimate", "Std. Error"))
    )
    expect_equal(coef(measlesFit, se = TRUE), orig,
                 tolerance = 1e-6) # increased for Solaris Sparc
    ## tolerance determined empirically by an R build with --disable-long-double
})

test_that("score vector and Fisher info agree with numerical approximations", {
    skip_if_not_installed("numDeriv")
    test <- function (neweights) {
        measlesModel$ne$weights <- neweights
        pencomp <- hhh4(measlesWeserEms, measlesModel,
                        check.analyticals = "numDeriv")$pen
        expect_equal(pencomp$score$analytic, pencomp$score$numeric,
                     tolerance = .Machine$double.eps^0.5)
        expect_equal(pencomp$fisher$analytic, pencomp$fisher$numeric,
                     tolerance = .Machine$double.eps^0.25)
    }
    test(W_powerlaw(maxlag = 5, normalize = FALSE, log = FALSE))
    ## normalized PL with maxlag < max(nbmat) failed in surveillance < 1.9.0:
    test(W_powerlaw(maxlag = 3, normalize = TRUE, log = TRUE))
})

test_that("automatic and manual normalization are equivalent", {
    ## check for equivalent functions
    for (type in c("powerlaw", "np")) {
        W_type <- get(paste0("W_", type), mode = "function")
        w0 <- W_type(maxlag = 3, normalize = TRUE)
        w1 <- surveillance:::scaleNEweights.list(
            W_type(maxlag = 3, normalize = FALSE),
            normalize = TRUE)
        pars <- w0$initial
        nbmat <- neighbourhood(measlesWeserEms)
        expect_equal(w1$w(pars, nbmat), w0$w(pars, nbmat))
        ## for the power law, dw and d2w are length 1 lists in w1 but not in w0
        unlistIfPL <- if (type == "powerlaw") function (x) x[[1L]] else identity
        expect_equal(unlistIfPL(w1$dw(pars, nbmat)), w0$dw(pars, nbmat))
        expect_equal(unlistIfPL(w1$d2w(pars, nbmat)), w0$d2w(pars, nbmat))
        ## microbenchmark::microbenchmark(w1$d2w(pars, nbmat), w0$d2w(pars, nbmat))
        ## -> type-specific implementations of normalized derivatives are faster
    }
    ## check for equivalent fits (rather redundant)
    measlesFit2 <- hhh4(
        stsObj = measlesWeserEms,
        control = modifyList(measlesModel, list(
            ne = list(
                weights = W_powerlaw(maxlag = 5, normalize = FALSE),
                normalize = TRUE # -> use scaleNEweights.list()
                )))
        )
    expect_equal(measlesFit, measlesFit2, ignore = "control",
                 tolerance = 1e-6) # increased to pass on 32-bit Windows
})
