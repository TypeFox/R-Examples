context("Unit test of the .xfcvl-familiy of functions")

lambda. <- symm(matrix(rchisq(9, df = pi), 3, 3))
ns. <- c(50, 60, 70)
Ylist. <- createS(n = ns., p = 15, dataset = TRUE, topology = "banded")
Tlist. <- replicate(length(ns.), rchisq(1,5)*diag(15), simplify = FALSE)
k. <- 10


for (i in 1:4) {

  res <- switch(i,
                rags2ridges:::.fcvl(lambda., Ylist., Tlist.),
                rags2ridges:::.sfcvl(lambda., Ylist., Tlist.),
                rags2ridges:::.kfcvl(lambda., Ylist., Tlist., k = k.),
                rags2ridges:::.afcvl(lambda., Ylist., Tlist.))

  test_that(".xfcl functions returns correctly formatted output", {
    expect_that(res, is_a("numeric"))
    expect_that(length(res), equals(1L))
  })

}

test_that(".xfcl functions works properly on degenerated data", {
  expect_that(TRUE, is_true())  # To be tested
})

# Expand tests







