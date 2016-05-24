# Test plotdelta function #

# Test that plotdelta outputs a ggplot object #
test_that("plotdelta produces a graph", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  testdata <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                         baseline = lm(Mass ~ 1, data = Mass), furthest = 3, closest = 2, 
                         type = "variable", stat = "max", func = "lin", cmissing = FALSE)
  
  
  test <- plotdelta(dataset = testdata[[1]]$Dataset)
  
  expect_true(attr(test, "class")[1] == "gg")
  
})