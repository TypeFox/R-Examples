# Test the plothist function #

# Test that plothist produces a ggplot object #
test_that("plothist produces a graph", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  testdata <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                         baseline = lm(Mass ~ 1, data = Mass), furthest = 3, closest = 2, 
                         type = "variable", stat = "max", func = "lin", cmissing = FALSE)
  
  testdatarand <- randwin(repeats = 2, xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                          baseline = lm(Mass ~ 1, data = Mass), furthest = 3, closest = 2, 
                          type = "variable", stat = "max", func = "lin", cmissing = FALSE)
  
  
  test  <- plothist(dataset = testdata[[1]]$Dataset)
  test2 <- plothist(dataset = testdata[[1]]$Dataset, datasetrand = testdatarand, histq = 0.95)
  
  expect_true(attr(test, "class")[1] == "gg")
  expect_true(attr(test2, "class")[1] == "gg")
  expect_false(is.null(test2$facet$facets))

})