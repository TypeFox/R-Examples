# Test function plotbetas #

# Test that plotbetas produces a ggplot objects #
# Test that plotbetas produces new plots for quad and cub #
test_that("plotbetas produces a graph", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  testdata <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                         baseline = lm(Mass ~ 1, data = Mass), furthest = 3, closest = 2, 
                         type = "variable", stat = "max", func = "lin", cmissing = FALSE)
  
  testenv <- environment()
  test    <- plotbetas(dataset = testdata[[1]]$Dataset)
  
  expect_true(attr(test, "class")[1] == "gg")
 
  testdata[[1]]$Dataset$ModelBetaQ <- testdata[[1]]$Dataset$ModelBeta
  testdata[[1]]$Dataset$Function   <- "quad"
  
  test <- plotbetas(dataset = testdata[[1]]$Dataset, plotall = TRUE, plotallenv = testenv)
  
  expect_true(exists("beta2", envir = testenv))
  
  testdata[[1]]$Dataset$ModelBetaC <- testdata[[1]]$Dataset$ModelBeta
  testdata[[1]]$Dataset$Function   <- "cub"
  
  test <- plotbetas(dataset = testdata[[1]]$Dataset, plotall = TRUE, plotallenv = testenv)
  
  expect_true(exists("beta2", envir = testenv))
  expect_true(exists("beta3", envir = testenv))
  
})