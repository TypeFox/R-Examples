# Test the PlotCor and PlotWeight functions #

# Test that PlotCor creates a ggplot object #
test_that("plotcor produces a graph", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  cross <- crosswin(xvar = list(Temp = MassClimate$Temp), 
                    xvar2 = list(Rain = MassClimate$Rain), cdate = MassClimate$Date,
                    bdate = Mass$Date, furthest = 2, closest = 1, 
                    stat = "max", stat2 = "max", type = "variable",
                    cmissing = FALSE, cinterval = "day")
  
  test <- plotcor(cor.output = cross, type = "A")
  
  expect_true(attr(test, "class")[1]=="gg") 
  
})