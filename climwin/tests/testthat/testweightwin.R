# Test weightwin function #

# Test that weightwin has output BestModel, BestModelData and WeightedOutput #
# Check that BestModel has coefficients not NAs
# Check there are no NAs in BestModelData or WeightedOutput
# Check that BestModelData is at least 2 columns
# Check that WeightedOutput has at least 7 columns
test_that("weightwin test", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- weightwin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date,
                    bdate = Mass$Date, baseline = lm(Mass ~ 1, data = Mass),
                    furthest = 2, closest = 1, func = "lin",
                    type = "variable", weightfun = "W", cinterval = "day",
                    par = c(3, 0.2, 0), control = list(ndeps = c(0.01, 0.01, 0.01)),
                    method = "L-BFGS-B")
  
  expect_error(weightwin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date,
                         bdate = Mass$Date, baseline = lm(Mass$Mass ~ 1),
                         furthest = 2, closest = 1, func = "lin",
                         type = "variable", weightfun = "W", cinterval = "day",
                         par = c(-1, 0.2, 0), control = list(ndeps = c(0.01, 0.01, 0.01)),
                         method = "L-BFGS-B"))
  
  expect_error(weightwin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date,
                         bdate = Mass$Date, baseline = lm(Mass$Mass ~ 1),
                         furthest = 2, closest = 1, func = "lin",
                         type = "variable", weightfun = "W", cinterval = "day",
                         par = c(3, -1, 0), control = list(ndeps = c(0.01, 0.01, 0.01)),
                         method = "L-BFGS-B"))
  
  expect_error(weightwin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date,
                         bdate = Mass$Date, baseline = lm(Mass$Mass ~ 1),
                         furthest = 2, closest = 1, func = "lin",
                         type = "variable", weightfun = "W", cinterval = "day",
                         par = c(3, 0.2, 1), control = list(ndeps = c(0.01, 0.01, 0.01)),
                         method = "L-BFGS-B"))
  
  expect_error(weightwin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date,
                         bdate = Mass$Date, baseline = lm(Mass$Mass ~ 1),
                         furthest = 2, closest = 1, func = "lin",
                         type = "variable", weightfun = "G", cinterval = "day",
                         par = c(3, -1, 0), control = list(ndeps = c(0.01, 0.01, 0.01)),
                         method = "L-BFGS-B"))
  

  expect_true(is.list(test))
  expect_false(is.na(test[[1]][1]))
  
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
  expect_true(is.list(test[[3]]))
  expect_equal(length(which(is.na(test[[3]]))), 0)
  
})
