## library(testthat)
context("pd0 and d.prime0 arguments to discrim")

##################################################################
## Tests of triangle family object:

test_that("Expect error if more than one of pd0/d.prime0 has been specified", {
    expect_error(
        discrim(26, 75, method = "triangle", d.prime0 = 2, pd0=.2, test =
                "similarity"),
        "Only specify one of")
    ## Expect error if none of pd0/d.prime0 has been specified:
    expect_error(
        discrim(26, 75, method = "triangle", test = "similarity")
        , "Either 'pd0' or 'd.prime0' has to be specified for a similarity test")
})

test_that("Specification of different scales give the same results:", {
T1 <- discrim(26, 75, method = "triangle", pd0 = .2, test = "similarity")
T2 <- discrim(26, 75, method = "triangle", test = "similarity",
              d.prime0 = psyinv(pd2pc(.2, 1/3), "triangle"))
expect_equal(T1$p.value, T2$p.value, tolerance=1e-3)
})

test_that("Test boundary values for d.prime0 (-1, 0, 1, Inf):", {
expect_that(
    discrim(26, 75, method = "triangle", d.prime0 = -1, test =
            "similarity")
    , throws_error("d.prime0 >= 0"))
expect_that(
    discrim(26, 75, method = "triangle", d.prime0 = 0, test =
            "similarity")
    , gives_warning("'d.prime0' should be positive for a similarity test"))

expect_output(
    print(discrim(26, 75, method = "triangle", d.prime0 = 1, test =
            "similarity"))
    , "p-value = 0.1274")

expect_output(
    print(discrim(26, 75, method = "triangle", d.prime0 = Inf,
            stat="like", test="similarity"))
    , "d-prime is less than Inf")
})

test_that("Test boundary values for pd0 (-1, 0, .2. 1. 2):", {
expect_that(
    discrim(26, 75, method = "triangle", pd0 = -1, test =
            "similarity")
    , throws_error("pd0 >= 0 is not TRUE"))
expect_that(
    discrim(26, 75, method = "triangle", pd0 = 0, test =
            "similarity")
    , gives_warning("'pd0' should be positive for a similarity test"))
expect_output(
    print(discrim(26, 75, method = "triangle", pd0 = .2, test =
            "similarity"))
    , "'exact' binomial test:  p-value = 0.02377")

expect_output(
    print(discrim(26, 75, method = "triangle", pd0 = 1, test = "similarity"))
    , "'exact' binomial test:  p-value = < 2.2e-16")

expect_error(
    discrim(26, 75, method = "triangle", pd0 = 2, test = "similarity")
    , "pd0 <= 1 is not TRUE")
})

test_that("Test that all statistics works in the limit of pd0 and d.prime0", {
Stat <- eval(formals(discrim)$statistic)
pvals <- sapply(Stat, function(stat) {
    discrim(26, 75, method = "triangle", d.prime0 = Inf,
            stat=stat, test="similarity")$p.value
})
expect_equivalent(pvals, rep(0, length(Stat)))

pvals <- sapply(Stat, function(stat) {
    discrim(26, 75, method = "triangle", d.prime0 = Inf,
            stat=stat, test="diff")$p.value
})
expect_equivalent(pvals, rep(1, length(Stat)))

pvals <- sapply(Stat, function(stat) {
    discrim(26, 75, method = "triangle", pd0=1,
            stat=stat, test="simi")$p.value
})
expect_equivalent(pvals, rep(0, length(Stat)))
pvals <- sapply(Stat, function(stat) {
    discrim(26, 75, method = "triangle", pd0=1,
            stat=stat, test="diff")$p.value
})
expect_equivalent(pvals, rep(1, length(Stat)))
})

test_that("Test error at invalid args for pd0 and d.prime0", {
expect_error( ## vector pd0
    discrim(26, 75, pd0=1:2)
    )
expect_error( ## character pd0
    discrim(26, 75, pd0="2")
    )
expect_error(
    discrim(26, 75, d.prime0=1:2)
    )
expect_error(
    discrim(26, 75, d.prime="2")
    )
expect_error( ## list pd0
    discrim(26, 75, pd0=list(1))
)
})

test_that("Printing alternative hypothesis in terms of d-prime (default):", {
expect_output(
    print(discrim(26, 75, method = "triangle"))
    , "Alternative hypothesis: d-prime is greater than 0 ")

expect_equal(
    discrim(26, 75, method = "triangle")$pd0
    , 0)
expect_equal(
    discrim(26, 75, method = "triangle")$alt.scale
    , "d-prime")
expect_equal(
    discrim(26, 75, method = "triangle", test="simil",
            pd0=.2)$alt.scale
    , "pd")
expect_equal(
    discrim(26, 75, method = "triangle", test="simil",
            d.prime0=2)$alt.scale
    , "d-prime")
})

test_that("Test that alternative hypothesis uses d.prime/pd", {
expect_output(
    print(discrim(26, 75, d.prime=0, method = "triangle"))
    , "Alternative hypothesis: d-prime is greater than 0 ")
expect_output(
    print(discrim(26, 75, pd0=0, method = "triangle"))
    , "Alternative hypothesis: pd is greater than 0 ")
})

