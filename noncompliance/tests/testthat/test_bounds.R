library(noncompliance)
context("ACE bounds and IV inequality checks")

test_that("IV inequality checks", {
  expect_equal(object = Check_IV_ineqs(158, 14, 0, 0, 52, 12, 23, 78),
               expected = TRUE)
  expect_equal(object = Check_IV_ineqs(c(158, 14, 0, 0, 52, 12, 23, 78)),
               expected = TRUE)
  expect_equal(object = as.numeric(Check_IV_ineqs(158, 14, 0, 0,
                                                  52, 12, 23, 78, TRUE)[[1]]),
               expected = c(0.9913, 0.3965, 0.4727, 0.1394),
               tolerance = 1e-6)
  expect_equal(object = Check_IV_ineqs(99, 1027, 30, 233, 84, 935, 31, 422),
               expected = TRUE)
  expect_equal(object = Check_IV_ineqs(c(99, 1027, 30, 233, 84, 935, 31, 422)),
               expected = TRUE)
  expect_equal(object = as.numeric(Check_IV_ineqs(99, 1027, 30, 233,
                                                  84, 935, 31, 422, TRUE)[[1]]),
               expected = c(0.7065, 0.7964, 0.3083, 0.1888),
               tolerance = 1e-6)
})

test_that("ACE bounds", {
  expect_equal(object = ACE_bounds(158, 14, 0, 0, 52, 12, 23, 78),
               expected = c(0.3913, 0.7792),
               tolerance = 1e-4)
  expect_equal(object = ACE_bounds(c(158, 14, 0, 0, 52, 12, 23, 78)),
               expected = c(0.3913, 0.7792),
               tolerance = 1e-4)
  expect_equal(object = ACE_bounds(99, 1027, 30, 233, 84, 935, 31, 422),
               expected = c(-0.6420, 0.2390),
               tolerance = 1e-4)
  expect_equal(object = ACE_bounds(c(99, 1027, 30, 233, 84, 935, 31, 422)),
               expected = c(-0.6420, 0.2390),
               tolerance = 1e-4)
})

test_that("ACDE bounds", {
  expect_equal(object = as.numeric(Check_ACDE_bounds(99, 1027, 30, 233,
                                                     84, 935, 31, 422)),
               expected = c(-0.0824, 0.0205, 0.0028, 0.1140),
               tolerance = 1e-3)
  expect_equal(object = as.numeric(Check_ACDE_bounds(c(99, 1027, 30, 233,
                                                       84, 935, 31, 422))),
               expected = c(-0.0824, 0.0205, 0.0028, 0.1140),
               tolerance = 1e-3)
  expect_equal(object = as.numeric(Check_ACDE_bounds(99, 1027, 30, 233,
                                                     84, 935, 31, 422,
                                                     iv.ineqs=TRUE)),
               expected = c(0.0142, 0.1042, 0.1189, -0.0005),
               tolerance = 1e-3)
  expect_equal(object = as.numeric(Check_ACDE_bounds(c(99, 1027, 30, 233,
                                                       84, 935, 31, 422),
                                                     iv.ineqs=TRUE)),
               expected = c(0.0142, 0.1042, 0.1189, -0.0005),
               tolerance = 1e-3)
})