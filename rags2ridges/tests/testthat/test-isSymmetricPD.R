context("Unit tests of the isSymmetricPD function")

# Alternative testing function
isSymmetricPDAlt <- function(M) {

  nm <- deparse(substitute(M))
  if (!is.matrix(M) || !is.numeric(M)) {
    stop(nm, " is not a numeric matrix")
  }
  if (!isSymmetric(M)) {
    stop(nm, " is not a symmetric matrix")
  }

  # Slower alternative test
  if (!all(eigen(M, symmetric = TRUE)$values > .Machine$double.eps)) {
    return(FALSE)
  } else {
    return(TRUE)
  }

}


test_that("isSymmetricPD works as intended", {

  pdS    <- createS(n = 15, p = 10)
  notpdS <- createS(n = 5, p = 10)

  expect_that(isSymmetricPD(pdS),    is_true())
  expect_that(isSymmetricPD(notpdS), is_false())

})


test_that("isSymmetricPD works for degenerate input", {

  # 0 by 0 matrices
  expect_false(isSymmetricPD(matrix(1,0,0)))

  # 1 by 1 matrices
  expect_true(isSymmetricPD(matrix(1,1,1)))
  expect_false(isSymmetricPD(matrix(0,1,1)))

})


test_that("isSymmetricPD throws errors when appropriate", {

  S1 <- createS(n = 15, p = 10)
  S2 <- createS(n = 5, p = 10)

  # Asymmetric
  S1[1,2] <- 1
  expect_that(isSymmetricPD(S1), throws_error("symmetric"))

  # Not a numeric matrix
  S2 <- as.character(S2)
  expect_that(isSymmetricPD(S2), throws_error("numeric"))

  # Not a matrix but numeric
  S2 <- as.numeric(S2)
  expect_that(isSymmetricPD(S2), throws_error("matrix"))

})


# test_that("isSymmetricPD agrees with isSymmtericPDAlt", {
#
#   for (n in round(seq(1, 20, l = 10))) {
#     for (p in round(seq(2, 20, l = 11))) {
#       cat(n, "    ",p, "\n")
#       S <- createS(n = n, p = p)
#       isSymmetricPD(S)
#       isSymmetricPDAlt(S)
#
#       expect_that(isSymmetricPD(S), equals(isSymmetricPDAlt(S)))
#     }
#   }
#
# })



