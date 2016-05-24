context("Neighbourhood order")

## generate random adjancency matrix
## radjmat <- function (n) {
##     adjmat <- matrix(0L, n, n, dimnames=list(letters[1:n],letters[1:n]))
##     adjmat[lower.tri(adjmat)] <- sample(0:1, n*(n-1)/2, replace=TRUE)
##     adjmat + t(adjmat)
## }
## set.seed(3); adjmat <- radjmat(5)
adjmat <- structure(
    c(0L, 0L, 1L, 0L, 0L,
      0L, 0L, 1L, 1L, 0L,
      1L, 1L, 0L, 0L, 1L,
      0L, 1L, 0L, 0L, 1L,
      0L, 0L, 1L, 1L, 0L),
    .Dim = c(5L, 5L),
    .Dimnames = rep.int(list(c("a", "b", "c", "d", "e")), 2L)
    )

## validated matrix of neighbourhood orders
nbmat <- structure(
    c(0L, 2L, 1L, 3L, 2L,
      2L, 0L, 1L, 1L, 2L,
      1L, 1L, 0L, 2L, 1L,
      3L, 1L, 2L, 0L, 1L,
      2L, 2L, 1L, 1L, 0L),
    .Dim = c(5L, 5L),
    .Dimnames = rep.int(list(c("a", "b", "c", "d", "e")), 2L)
    )

test_that("nbOrder() returns the validated matrix", {
    skip_if_not_installed("spdep")
    expect_that(suppressMessages(nbOrder(adjmat, maxlag=Inf)),
                is_identical_to(nbmat))
})

test_that("zetaweights(.,maxlag=1,normalize=FALSE) is inverse of nbOrder", {
    expect_that(zetaweights(nbmat, maxlag=1, normalize=FALSE),
                is_identical_to(1*adjmat))
})
