require(testthat)
require(pez)
require(caper)
data(laja)
data <- comparative.comm(invert.tree, river.sites, invert.traits, warn=FALSE)

context("phy.signal")

test_that("lambda", {
    set.seed(123)
    lambda <- phy.signal(data, method="lambda")
    expect_that(names(lambda), equals(names(data$data)))
    expect_that(round(lambda,3), is_equivalent_to(c(1, 0.294)))
})

test_that("delta", {
    set.seed(123)
    delta <- phy.signal(data, method="delta")
    expect_that(names(delta), equals(names(data$data)))
    expect_that(round(delta,3), is_equivalent_to(c(0.387, 3)))
})

test_that("kappa", {
    set.seed(123)
    kappa <- phy.signal(data, method="kappa")
    expect_that(names(kappa), equals(names(data$data)))
    expect_that(round(kappa,3), is_equivalent_to(c(0,0.061)))
})
