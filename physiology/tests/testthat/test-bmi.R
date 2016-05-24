context("ideal weight")

#TODO: common tests for all functions with invalid, severely out-of-range inputs

test_that("ideal_weight_adultAdult", {
  inch <- 100 / 2.54
  expect_error(ideal_weight_adult(male = TRUE))
  expect_error(ideal_weight_adult(heightm = 1.7))

  # should warn when height is out of validated range of the formula
  expect_warning(ideal_weight_adult(heightm = 0, male = TRUE, warn = TRUE))
  # should warn when height is out of validated range of the formula
  expect_warning(ideal_weight_adult(heightm = -1, male = TRUE, warn = TRUE))
  # should warn when height is out of validated range of the formula
  expect_warning(ideal_weight_adult(heightm = 3, male = TRUE, warn = TRUE))
  # should warn when height is out of validated range of the formula
  expect_warning(ideal_weight_adult(heightm = 59 / inch,
                                    male = TRUE,
                                    warn = TRUE))
  expect_that(ideal_weight_adult(heightm = 59 / inch,
                                    male = TRUE,
                                    warn = FALSE),
              testthat::not(gives_warning()))


  expect_equal(ideal_weight_adult(60 / inch, male = TRUE), 50)
  expect_equal(ideal_weight_adult(60 / inch, male = F), 45.5)
  expect_equal(ideal_weight_adult(c(60 / inch, 60 / inch), male=c(FALSE, TRUE)),
               c(45.5, 50))
  expect_equal(ideal_weight_adult(c(60 / inch, 60 / inch, NA),
                                male=c(FALSE, TRUE, TRUE)),
               c(45.5, 50, NA))
  expect_equal(ideal_weight_adult(c(60 / inch, 60 / inch, 60 / inch),
                                male = c(FALSE, NA, TRUE)),
               c(45.5, NA, 50))
  expect_error(ideal_weight_adult(c(60 / inch, 60 / inch, 60 / inch),
                                male = c(FALSE, TRUE)))
  expect_error(ideal_weight_adult(c(60 / inch, 60 / inch),
                                male = c(FALSE, TRUE, TRUE)))
  expect_error(ideal_weight_adult(c(), male=c(FALSE, TRUE, TRUE)))
  expect_error(ideal_weight_adult(c(60 / inch, 60 / inch), male = c()))

  expect_warning(ideal_weight_adult(12 * 8.1 / inch, male = TRUE, warn = TRUE))

})
