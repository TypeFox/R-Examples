# Test randwin function #

# Test that randwin create only WindowOutputRand, not BestModel or BestModelData
# Test that WindowOutputRand is at least 13 columns (will vary with func = Q or C)
# Test that the number of rows in each repeat is the same as number of windows
# Test that the randimised column is yes

test_that("Check randwin output", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  furthest = 2
  closest = 1
  stat = "max"
  
  rand <- randwin(repeats = 1, xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date, 
                  bdate = Mass$Date, baseline = lm(Mass ~ 1, data = Mass), 
                  furthest = 2, closest = 1, 
                  type = "variable", stat = "max", func = "lin", cmissing = FALSE)
  
  duration  <- (furthest - closest) + 1
  maxmodno  <- (duration * (duration + 1)) / 2
  
  expect_true(is.data.frame(rand))
  expect_equal(length(which(is.na(rand[, 4]))), 0)
  expect_true(ncol(rand) >= 15)
  expect_equal(maxmodno, nrow(subset(rand, Repeat == 1)))
  expect_true((rand["Randomised"])[1,] == "yes")
  
})