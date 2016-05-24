# Test the outcome of crosswin #

# Test that crosswin has created a correct CrossWinOutput object
# Expect that an object CrossWinOutput exists
# Expect that there are no NA values
# Expect that the number of columns is at least 7 (will vary with values of FIXED) 
# Expect that the number of rows is equal to the number of possible windows
test_that("crosswin produces output", {
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- crosswin(xvar = list(Temp = MassClimate$Temp), 
                   xvar2 = list(Rain = MassClimate$Rain), 
                   cdate = MassClimate$Date,
                   bdate = Mass$Date, furthest = 2, closest = 1, 
                   stat = "max", stat2 = "max", type = "variable",
                   cmissing = FALSE, cinterval = "day")
  
  furthest = 2
  closest = 1
  duration  <- (furthest - closest) + 1
  maxmodno  <- (duration * (duration + 1))/2
  
  expect_true(is.data.frame(test))
  expect_equal(length(which(is.na(test))), 0)
  expect_true(ncol(test) >= 7)
  expect_equal(maxmodno, nrow(test))
  
})