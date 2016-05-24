context("Checking the behavior of small utility functions")

test_that("Euler rotation returns a matrix", {
  m <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  m2 <- rbind(c(0, 1, 0), c(-1, 0, 0), c(0, 0, 1))
  expect_equal(cda$euler(0, 0, 0), m)
  expect_equal(cda$euler(pi/2, 0, 0), m2)
})

test_that("Dielectric function is correct", {
  .gold <- structure(list(wavelength = c(400, 500, 600), 
                          epsilon = c(-1.64965688407203+5.77176308089817i, 
                                      -2.99223400282977+3.63036432730633i, 
                                      -9.07165375859948+1.40807251911031i)), 
                     .Names = c("wavelength", "epsilon"), row.names = c(NA, -3L), 
                     class = "data.frame")
  gold <- epsAu(seq(400, 600, by=100))
  expect_equal(gold, .gold)
})
