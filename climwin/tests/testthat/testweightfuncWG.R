# Test loglik functions G and W #

# Test GEV function #
# Test that deltaAICc output is not NA #
# Test that deltaAICc is less than 0 #
test_that("GEV loglik test", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  funcenv              <- environment()
  xvar                 <- MassClimate$Temp
  bdate                <- Mass$Date
  cdate                <- MassClimate$Date
  furthest             <- 2
  closest              <- 1
  duration             <- (furthest-closest) + 1
  baseline             <- lm(Mass ~ 1, data = Mass)
  cmatrix              <- matrix(ncol = (duration), nrow = length(bdate))
  nullmodel            <- AICc(baseline)
  funcenv$modno        <- 1
  funcenv$DAICc        <- list()
  funcenv$par_shape    <- list()
  funcenv$par_scale    <- list()
  funcenv$par_location <- list()
  
  cont <- convertdate(bdate = bdate, cdate = cdate, xvar = xvar, 
                      cinterval = "day", type = "variable")   # create new climate dataframe with continuous daynumbers, leap days are not a problem 
  
  for (i in 1:length(bdate)){
    for (j in closest:furthest){
      k <- j - closest + 1
      cmatrix[i, k] <- xvar[match(cont$bintno[i] - j,cont$cintno)]   #Create a matrix which contains the climate data from furthest to furthest from each biological record#    
    }
  }
  
  modeldat         <- model.frame(baseline)
  modeldat$climate <- matrix(ncol = 1, nrow = nrow(cmatrix), seq(from = 1, to = nrow(cmatrix), by = 1))
  
  test <- modloglik_G(par = c(3, 0.2, 0), 
                      modeloutput = lm(Mass ~ climate, data = modeldat), 
                      duration = duration, cmatrix = cmatrix, 
                      nullmodel = nullmodel, funcenv = funcenv)
  
  expect_false(is.na(test))
  expect_true(test <= 2)
  
})

# Test Weibull function #
# Test that deltaAICc output is not NA #
# Test that deltaAICc is less than 0 #
test_that("Weibull loglik test", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  funcenv              <- environment()
  xvar                 <- MassClimate$Temp
  bdate                <- Mass$Date
  cdate                <- MassClimate$Date
  furthest             <- 2
  closest              <- 1
  duration             <- (furthest-closest) + 1
  baseline             <- lm(Mass ~ 1, data = Mass)
  cmatrix              <- matrix(ncol = (duration), nrow = length(bdate))
  nullmodel            <- AICc(baseline)
  funcenv$modno        <- 1
  funcenv$DAICc        <- list()
  funcenv$par_shape    <- list()
  funcenv$par_scale    <- list()
  funcenv$par_location <- list()
  
  cont <- convertdate(bdate = bdate, cdate = cdate, xvar = xvar, 
                      cinterval = "day", type = "variable")   # create new climate dataframe with continuous daynumbers, leap days are not a problem 
  
  for (i in 1:length(bdate)){
    for (j in closest:furthest){
      k <- j - closest + 1
      cmatrix[i, k] <- xvar[match(cont$bintno[i] - j,cont$cintno)]   #Create a matrix which contains the climate data from furthest to furthest from each biological record#    
    }
  }
  
  modeldat         <- model.frame(baseline)
  modeldat$climate <- matrix(ncol = 1, nrow = nrow(cmatrix), seq(from = 1, to = nrow(cmatrix), by = 1))
  
  test <- modloglik_W(par = c(3, 0.2, 0), 
                      modeloutput = lm(Mass ~ climate, data = modeldat), 
                      duration = duration, cmatrix = cmatrix, 
                      nullmodel = nullmodel, funcenv = funcenv)
  
  expect_false(is.na(test))
  expect_true(test <= 2)
  
})