testthat::context("Test OLC shortening and recovery")

# https://github.com/google/open-location-code/blob/master/test_data/shortCodeTests.csv
testthat::test_that("OLC shortening works for simple cases", {
  testthat::expect_equal(shorten_olc("9C3W9QCJ+2VX",51.3701125,-1.217765625), "+2VX")
})

testthat::test_that("OLC shortening works +/- .000755", {
  testthat::expect_equal(shorten_olc("9C3W9QCJ+2VX",51.3708675,-1.217765625), "CJ+2VX")
  testthat::expect_equal(shorten_olc("9C3W9QCJ+2VX",51.3701125,-1.217010625), "CJ+2VX")
})

testthat::test_that("OLC shortening works +/- .0151", {
  testthat::expect_equal(shorten_olc("9C3W9QCJ+2VX",51.3852125,-1.217765625), "9QCJ+2VX")
  testthat::expect_equal(shorten_olc("9C3W9QCJ+2VX",51.3550125,-1.217765625), "9QCJ+2VX")
  testthat::expect_equal(shorten_olc("9C3W9QCJ+2VX",51.3701125,-1.232865625), "9QCJ+2VX")
  testthat::expect_equal(shorten_olc("9C3W9QCJ+2VX",51.3701125,-1.202665625), "9QCJ+2VX")
})

testthat::test_that("OLC recovery works",{
  testthat::expect_equal(recover_olc("9G8F+6X", 47.4, 8.6), "8FVC9G8F+6X")
})

testthat::test_that("The floating-point bug in the reference implementation does not regress", {
  testthat::expect_equal(recover_olc("22+", 42.899, 9.012), "8FJFW222+")
  testthat::expect_equal(recover_olc("22+", 14.95125, -23.5001), "796RXG22+")
})

testthat::test_that("NA support works", {
  testthat::expect_true(is.na(recover_olc("22+", 14.95125, NA)))
  testthat::expect_true(is.na(recover_olc("22+", NA, 14.95125)))
  testthat::expect_true(is.na(recover_olc(NA, 14.95125, 14.95125)))

})
