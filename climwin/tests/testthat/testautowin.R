context("autowin function")

# Test the outcomes of autowin #

# Test that autowin has created a correct AutoWinOutput object
# Expect that an object AutoWinOutput exists
# Expect that there are no NA values
# Expect that the number of columns is at least 7 (will vary with values of FIXED) 
# Expect that the number of rows is equal to the number of possible windows
test_that("AutoWinOutput has created an output", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  single <- singlewin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date,
                      baseline = lm(Mass ~ 1, data = Mass), furthest = 1, closest = 1,
                      stat = "mean", func = "lin",
                      type = "variable", cmissing = FALSE, cinterval = "day")
  
  test <- autowin(reference = single$BestModelData$climate,
                  xvar  = list(Temp = MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date,
                  baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 1, 
                  stat = "mean", func = "lin", type = "variable", cmissing = FALSE, cinterval = "day")
  
  
  furthest <- 2
  closest  <- 1
  duration <- (furthest - closest) + 1
  maxmodno <- (duration * (duration + 1))/2
  
  expect_true(exists("test"))
  expect_equal(length(which(is.na(test))), 0)
  expect_true(ncol(test) >= 7)
  expect_equal(maxmodno, nrow(test))
})