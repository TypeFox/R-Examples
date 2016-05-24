context("Unit test of armaRidgeP")

# To make the R versions of armaRidgePAnyTarget and armaRidgePScalarTarget
# available.
example("armaRidgeP", package = "rags2ridges",
        character.only = TRUE, echo = FALSE)

# The functions to test
armaRidgeP <- rags2ridges:::.armaRidgeP  # To avoid writing rags2ridges:::
aRidgePAnyTarget <- rags2ridges:::.armaRidgePAnyTarget
aRidgePScalarTarget <- rags2ridges:::.armaRidgePScalarTarget

# Values to test
test.lambdas <- c(1e-200, 1e-100, 1e-50, 1e-14, 1e-10, 1,
                  1e10, 1e50, 1e100, 1e200, 1e300, 1e500, Inf)
tgt.types <- c("DAIE", "DIAES", "DUPV", "DAPV", "DCPV", "DEPV", "Null")

#
# Test that the C++ version agree with the R implementations in
# help("armaRidgeP")
#

S <- unname(createS(n = 5, p = 10)) # Create some data
for (type in tgt.types) {
  tgt <- default.target(S, type = type, const = 1)
  for (j in 1:2) {
    a <- switch(j, "aRidgePAnyTarget", "aRidgePScalarTarget")
    r <- switch(j, "rRidgePAnyTarget", "rRidgePScalarTarget")
    t <- switch(j, tgt, tgt[1,1])
    for (l in c(1e-14, 1e-5, 1, 10, 1e4)) {
      for (invert in 0:2) {
        test_that(sprintf("%s() agrees with %s() for l=%g, type=%s, invert=%d",
                          a, r, l, type, invert), {
          expect_equal(get(a)(S, t, l, invert), get(r)(S, t, l, invert))

        })
      }
    }
  }
}

#
# Futher tests of armaRidgeP
#

p <- 4
n <- 5
for (n in c(5, 9, 14)) {
for (p in c(4, 10, 15)){

# Create some toy data
S <- unname(createS(n = n, p = p))

for (type in tgt.types) {
  tgt <- default.target(S, type = type, const = 1)

  for (l in test.lambdas) {

    if (type == "DEPV" && l <= 1e-50) {
      next
    }

    res <- armaRidgeP(S, tgt, l)

    test_that(paste("proper format for lambda =", l), {
      expect_that(is.double(res), is_true())  # Returns numeric (dobule)
      expect_that(res, is_a("matrix"))        # Returns a matrix
      expect_that(dim(res), equals(dim(S)))   # .. of the correct size
    })

  } ## End for l

  test_that(paste("proper values for very large lambda, tgt =", type), {

    expect_that(armaRidgeP(S, tgt, 1e200), equals(tgt))
    expect_that(armaRidgeP(S, tgt, Inf), equals(tgt))

  })

  test_that(paste("proper values for very small lambda, type =", type), {

    expect_that(armaRidgeP(S, tgt, 1e-10), not(equals(tgt)))
    expect_that(armaRidgeP(S, tgt, 1e-50), not(equals(tgt)))
    expect_that(armaRidgeP(S, tgt, 1e-100), not(equals(tgt)))
    expect_that(armaRidgeP(S, tgt, 1e-200), not(equals(tgt)))
    expect_that(armaRidgeP(S, tgt, 1e-300), not(equals(tgt)))
    expect_that(armaRidgeP(S, tgt, 1e-400), throws_error("postive"))
    expect_that(armaRidgeP(S, tgt, 0),      throws_error("postive"))

    if (p > n) {
      aa <- armaRidgeP(S, tgt, 1e-10)
      bb <- armaRidgeP(S, tgt, 1e-50)
      cc <- armaRidgeP(S, tgt, 1e-100)
      dd <- armaRidgeP(S, tgt, 1e-200)
      ee <- armaRidgeP(S, tgt, 1e-300)

      expect_that(all(abs(aa) <= abs(bb)), is_true())
      expect_that(all(abs(bb) <= abs(cc)), is_true())
      expect_that(all(abs(cc) <= abs(dd)), is_true())
      expect_that(all(abs(dd) <= abs(ee)), is_true())
    }

  })

} ## End for type

} ## End for p
} ## End for n


#
# Test for very large lambda AND targets
#

# source("../tests/testthat/reference-values.R")
source("reference-values.R")

test_that("Test armaRidgeP in various special cases (by reference)", {

  expect_that(any(!is.finite(armaRidgeP(Sbar, Tbar, aa))), is_false())
  expect_that(armaRidgeP(Sbar, Tbar, aa), equals(Tbar))

})
