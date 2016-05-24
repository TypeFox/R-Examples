# Test function plotbest #

# Test that plotbest produces a ggplot object #
test_that("plotbest produces a graph", {
  
  data(MassOutput, envir = environment())
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  testdata <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                         baseline = lm(Mass ~ 1, data = Mass), furthest = 3, closest = 2, 
                         type = "variable", stat = "max", func = "lin", cmissing = FALSE)
  
  single <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                      cdate = MassClimate$Date, bdate = Mass$Date,
                      baseline = lm(Mass$Mass ~ 1),furthest = 72, closest = 15,
                      stat = "max", func = "lin",
                      type = "variable", cmissing = FALSE, cinterval = "day")
  
  test <- plotbest(dataset = testdata[[1]]$Dataset, bestmodel = single[[1]],
           bestmodeldata = single[[2]])
  
  expect_true(attr(test, "class")[1] == "gg")
  
})