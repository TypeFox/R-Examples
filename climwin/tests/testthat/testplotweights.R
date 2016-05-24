# Test the plotweights function #

# Test that plotweights produces a ggplot object #
test_that("plotweights produces a graph", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  testdata <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                         baseline = lm(Mass ~ 1, data = Mass), furthest = 3, closest = 2, 
                         type = "variable", stat = "max", func = "lin", cmissing = FALSE)
  
  
  test <- plotweights(dataset = testdata[[1]]$Dataset, cw1 = 0.95, cw2 = 0.75, cw3 = 0.25)
  
  expect_true(attr(test, "class")[1] == "gg")
  
})