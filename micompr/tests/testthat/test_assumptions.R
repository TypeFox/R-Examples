# Copyright (c) 2016 Nuno Fachada
# Distributed under the MIT License (http://opensource.org/licenses/MIT)

library(micompr)
context("assumptions")

test_that("assumptions_manova constructs the expected objects", {


  # Test depends on whether the MVN and biotools packages are present
  if (!requireNamespace("MVN", quietly = TRUE) ||
      !requireNamespace("biotools", quietly = TRUE) ) {

    # MVN and biotools packages are NOT present

    expect_message(assumptions_manova(iris[, 1:4], iris[, 5]),
                   "MANOVA assumptions require 'MVN' and 'biotools' packages.")

    expect_null(assumptions_manova(iris[, 1:4], iris[, 5]))

  } else {

    # MVN and biotools packages are present

    amnv <- assumptions_manova(iris[, 1:4], iris[, 5])

    expect_is(amnv, "assumptions_manova")
    for (rt in amnv$mvntest) {
      expect_is(rt, "royston")
    }
    expect_is(amnv$vartest, "boxM")

  }

})

test_that("assumptions_manova throws the expected warnings", {

  # Create bogus data for testing
  bogus_data <- matrix(rnorm(100), ncol = 10)
  factors <- c(rep("A", 3), rep("B", 3), rep("C", 4))

  # Test depends on whether the MVN and biotools packages are present
  if (!requireNamespace("MVN", quietly = TRUE) ||
      !requireNamespace("biotools", quietly = TRUE) ) {

    # MVN and biotools packages are NOT present

    expect_message(assumptions_manova(bogus_data, factors),
                   "MANOVA assumptions require 'MVN' and 'biotools' packages.")

    expect_null(assumptions_manova(bogus_data, factors))

  } else {

    # MVN and biotools packages are present, perform test

    # This warning is thrown twice because of groups A and B
    expect_warning(assumptions_manova(bogus_data, factors),
                   "Royston test requires at least 4 observations",
                   fixed = TRUE)

    # This warning is thrown because of group C, which has more variables (10)
    # than observations (4)
    expect_warning(assumptions_manova(bogus_data, factors),
                   "Royston test requires more observations than (dependent)",
                   fixed = TRUE)


  }

})

test_that("assumptions_paruv constructs the expected objects", {

  auv <- assumptions_paruv(iris[, 1:4], iris[, 5])

  expect_is(auv, "assumptions_paruv")

  # Test that the Shapiro-Wilk test is present for all tested dependent
  # variables and groups
  for (dv in auv$uvntest) {

    for (swt in dv) {

      expect_is(swt, "htest")
      expect_output(print(swt),
                    "Shapiro-Wilk normality test",
                    fixed = TRUE)
    }
  }

  # Test that the Bartlett test is present for all tested dependent variables
  for (bt in auv$vartest) {

    expect_is(bt, "htest")
    expect_output(print(bt),
                  "Bartlett test of homogeneity of variances",
                  fixed = TRUE)
  }

})