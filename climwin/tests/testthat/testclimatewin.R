# Test the outcomes of climatewin #

#Need to vary:
# baseline model type: glm, lme4, glmer
# cinterval: D, W, M
# stat: mean, min, max, slope
# func: L, Q, C, I, LOG
# upper: range
# lower: range
# thresh: TRUE or FALSE
# cvk: 0 and >1
# centre: value

# Test that climatewin has created a BestModel, BestModelData and WindowOutput
# Expect that an object BestModel exists
# Expect that the coefficients of BestModel are not NAs
# Check there are no NAs in BestModelData or WindowOutput
# Check there are at least 2 columns in BestModelData
# Check there are at least 15 columns in WindowOutput
# Check that the number of rows in the same as the number of windows
# Check that randomised is "no"

##########################################################

# Test regular output

test_that("climatewin produces the right output", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  furthest = 2
  closest = 2
  
  test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                     type = "variable", stat = "max", func = "lin", cmissing = FALSE)
  
  duration  <- (furthest - closest) + 1
  maxmodno  <- (duration * (duration + 1))/2
  
  expect_true(is.list(test))
  expect_false(is.na((test[[1]]$BestModel)[1]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 2)
  
  expect_equal(length(which(is.na(test[[1]]$Dataset[, 4]))), 0)
  expect_true(ncol(test[[1]]$Dataset) >= 15)
  expect_equal(maxmodno, nrow(test[[1]]$Dataset))
  expect_true((test[[1]]$Dataset["Randomised"])[1, ] == "no")
  
})

test_that("climatewin produces multiple combos", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  furthest = 1
  closest = 0
  
  test <- climatewin(xvar = list(Temp = MassClimate$Temp, Rain = MassClimate$Rain), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 1, closest = 0, 
                     type = "variable", stat = c("max", "min"), func = c("lin", "quad"), cmissing = FALSE)
  
  expect_equal(nrow(test$combos), 8)
  
})

##########################################################

# Test that upper and lower work

test_that("climatewin produces binary values with upper and thresh = TRUE", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  furthest = 1
  closest = 0
  
  test <- climatewin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 1, closest = 0, 
                     type = "variable", stat = "max", func = "lin", cmissing = FALSE,
                     upper = 10, thresh = TRUE)
  
  expect_equal(min(test[[1]]$BestModelData$climate), 0)
  expect_equal(max(test[[1]]$BestModelData$climate), 1)
  
})

test_that("climatewin produces non-binary values with upper and thresh = FALSE", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  furthest = 1
  closest = 0
  
  test <- climatewin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 1, closest = 0, 
                     type = "variable", stat = "max", func = "lin", cmissing = FALSE,
                     upper = 10, thresh = FALSE)
  
  expect_equal(min(test[[1]]$BestModelData$climate), 0)
  expect_true(max(test[[1]]$BestModelData$climate) > 1)
  
})

test_that("climatewin produces non-binary values with lower and thresh = FALSE", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  furthest = 1
  closest = 0
  
  test <- climatewin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 1, closest = 0, 
                     type = "variable", stat = "max", func = "lin", cmissing = FALSE,
                     lower = 10, upper = 15, thresh = FALSE)
  
  expect_equal(min(test[[1]]$BestModelData$climate), 0)
  expect_true(max(test[[1]]$BestModelData$climate) > 1)
  
})

test_that("climatewin produces binary values with lower and thresh = TRUE", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  furthest = 1
  closest = 0
  
  test <- climatewin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 1, closest = 0, 
                     type = "variable", stat = "max", func = "lin", cmissing = FALSE,
                     lower = 10, thresh = TRUE)
  
  expect_equal(min(test[[1]]$BestModelData$climate), 0)
  expect_equal(max(test[[1]]$BestModelData$climate), 1)
  
})

test_that("climatewin produces binary values with lower/upper and thresh = TRUE", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  furthest = 1
  closest = 0
  
  test <- climatewin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 1, closest = 0, 
                     type = "variable", stat = "max", func = "lin", cmissing = FALSE,
                     lower = 10, upper = 15, thresh = TRUE)
  
  expect_equal(min(test[[1]]$BestModelData$climate), 0)
  expect_equal(max(test[[1]]$BestModelData$climate), 1)
  
})

test_that("climatewin produces non-binary values with lower/upper and thresh = FALSE", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  furthest = 1
  closest = 0
  
  test <- climatewin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 1, closest = 0, 
                     type = "variable", stat = "max", func = "lin", cmissing = FALSE,
                     lower = 10, upper = 15, thresh = FALSE)
  
  expect_equal(min(test[[1]]$BestModelData$climate), 0)
  expect_true(max(test[[1]]$BestModelData$climate) > 1)
  
})

##########################################################

# Test different settings of cmissing #

#When cmissing is TRUE and no NA is present#
test_that("No errors return when cmissing TRUE and full dataset", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                     type = "variable", stat = "max", func = "lin", cmissing=TRUE)
  
  expect_true(is.list(test))

})

#When cmissing is TRUE and NA is present#
test_that("No errors return when cmissing TRUE with NAs", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  MassClimate2 <- MassClimate[-491, ]
  test <- climatewin(xvar = list(MassClimate2$Temp), cdate = MassClimate2$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 0, 
                     type = "variable", stat = "max", func = "lin", cmissing = TRUE)
  
  expect_true(is.list(test))
  
})

#When cmissing is FALSE and NA is present#
#Test that an error occurs
#Test that object missing is made
#Test that object missing has length 1 (only 1 Date has been removed)

test_that("Error returned when cmissing FALSE with NAs, cinterval = D", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  MassClimate2 <- MassClimate[-491, ]
  expect_error(climatewin(xvar = list(MassClimate2$Temp), cdate = MassClimate2$Date, bdate = Mass$Date, 
                          baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                          type = "variable", stat = "max", func = "lin", 
                          cmissing=FALSE))
  
  expect_true(exists("missing"))
  expect_equal(length(missing), 1)
  rm("missing", envir = .GlobalEnv)
  
  })

test_that("Error returned when cmissing FALSE with NAs, cinterval = W", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  MassClimate2 <- MassClimate[-491, ]
  expect_error(climatewin(xvar = list(MassClimate2$Temp), cdate = MassClimate2$Date, bdate = Mass$Date, 
                          baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 0, 
                          type = "variable", stat = "max", func = "lin", cinterval = "week",
                          cmissing=FALSE))
  
  expect_true(exists("missing"))
  expect_equal(length(missing), 1)
  rm("missing", envir = .GlobalEnv)
  
})

test_that("Error returned when cmissing FALSE with NAs, cinterval = M", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  MassClimate2 <- MassClimate[-491, ]
  expect_error(climatewin(xvar = list(MassClimate2$Temp), cdate = MassClimate2$Date, bdate = Mass$Date, 
                          baseline = lm(Mass ~ 1, data = Mass), furthest = 1, closest = 0, 
                          type = "variable", stat = "max", func = "lin", cinterval = "month",
                          cmissing=FALSE))
  
  expect_true(exists("missing"))
  expect_equal(length(missing), 1)
  rm("missing", envir = .GlobalEnv)
  
})

##########################################################

# Test different types of models #

# Test glm models #
test_that("glm models can run", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = glm(Mass ~ 1, data = Mass, family = poisson), furthest = 2, closest = 2, 
                     type = "variable", stat = "max", func = "lin", cmissing=FALSE)
  
  expect_true(is.list(test))
  expect_false(is.na((test[[1]]$BestModel)[1]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 2)
  
})

test_that("lmer models can run", {
  
  data(Offspring, envir = environment())
  data(OffspringClimate, envir = environment())
  
  test <- climatewin(xvar = list(OffspringClimate$Temp), cdate = OffspringClimate$Date, 
                     bdate = Offspring$Date, 
                     baseline = lmer(Offspring ~ 1 + (1|BirdID), data = Offspring),  
                     furthest = 2, closest = 2, type = "variable", 
                     stat = "max", func = "lin", cmissing=FALSE)
  
  expect_true(is.list(test))
  expect_false(is.na(fixef(test[[1]]$BestModel)[1]))
  expect_false(is.na(fixef(test[[1]]$BestModel)[2]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 2)
  
})

test_that("glmer models can run", {
  
  data(Offspring, envir = environment())
  data(OffspringClimate, envir = environment())
  
  suppressWarnings(test <- climatewin(xvar = list(OffspringClimate$Temp), cdate = OffspringClimate$Date, 
                     bdate = Offspring$Date, 
                     baseline = glmer(Offspring ~ 1 + (1|Order), data = Offspring, family = "poisson"),  
                     furthest = 1, closest = 0, type = "variable", 
                     stat = "max", func = "lin", cmissing=FALSE))
  
  expect_true(is.list(test))
  expect_false(is.na(fixef(test[[1]]$BestModel)[1]))
  expect_false(is.na(fixef(test[[1]]$BestModel)[2]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 2)
  
})

##########################################################

# Test fixed and variable windows#

# Test fixed window#
test_that("Fixed window works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                     type = "fixed", cutoff.day = 20, cutoff.month = 5, 
                     stat = "max", func = "lin", cmissing=FALSE)
  
  expect_true(is.list(test))
  expect_false(is.na((test[[1]]$BestModel)[1]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 2)
  
})

##########################################################

# Test slope stat #
test_that("slope stat work", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 1, 
                     type = "variable", stat = "slope", func = "lin", cmissing=FALSE)

  expect_true(is.list(test))
  expect_false(is.na((test[[1]]$BestModel)[1]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 2)
  
})

test_that("slope and LOG return error", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  expect_error(cliamtewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date,
                          bdate = Mass$Date, baseline = lm(Mass ~ 1, data = Mass),
                          furthest = 2, closest = 1, type = "variable", stat = "slope",
                          func = "log"))  
  
})

test_that("slope and I return error", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  expect_error(cliamtewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date,
                          bdate = Mass$Date, baseline = lm(Mass ~ 1, data = Mass),
                          furthest = 2, closest = 1, type = "variable", stat = "slope",
                          func = "inv"))  
  
})

##########################################################

#Test different functions for fitting climate#

#Test quadratic function#
test_that("Quadratic function works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                     type = "variable", stat = "max", func = "quad", cmissing=FALSE)
  
  furthest = 2
  closest = 2
  duration  <- (furthest - closest) + 1
  maxmodno  <- (duration * (duration + 1))/2
  
  expect_true(is.list(test))
  expect_false(is.na((test[[1]]$BestModel)[1]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 3)
  
  expect_equal(length(which(is.na(test[[1]]$Dataset[, 5]))), 0)
  expect_true(test[[1]]$Dataset[, 8] == "quad")
  expect_true(ncol(test[[1]]$Dataset) >= 15)
  expect_equal(maxmodno, nrow(test[[1]]$Dataset))
  expect_true((test[[1]]$Dataset["Randomised"])[1, ] == "no")
  
})

#Test cubic function#
test_that("Cubic function works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                     type = "variable", stat = "max", func = "cub", cmissing=FALSE)
  
  furthest = 2
  closest = 2
  duration  <- (furthest - closest) + 1
  maxmodno  <- (duration * (duration + 1))/2
  
  expect_true(is.list(test))
  expect_false(is.na((test[[1]]$BestModel)[1]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 4)
  
  expect_equal(length(which(is.na(test[[1]]$Dataset[, 6]))), 0)
  expect_true(test[[1]]$Dataset[, 8] == "cub")
  expect_true(ncol(test[[1]]$Dataset) >= 15)
  expect_equal(maxmodno, nrow(test[[1]]$Dataset))
  expect_true((test[[1]]$Dataset["Randomised"])[1, ] == "no")
 
})

#Test log function#
test_that("Log function works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                     type = "variable", stat = "max", func = "log", cmissing=FALSE)
  
  furthest = 2
  closest = 2
  duration  <- (furthest - closest) + 1
  maxmodno  <- (duration * (duration + 1))/2
  
  expect_true(is.list(test))
  expect_false(is.na((test[[1]]$BestModel)[1]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 2)
  
  expect_equal(length(which(is.na(test[[1]]$Dataset[, 4]))), 0)
  expect_true(test[[1]]$Dataset[, 8] == "log")
  expect_true(ncol(test[[1]]$Dataset) >= 15)
  expect_equal(maxmodno, nrow(test[[1]]$Dataset))
  expect_true((test[[1]]$Dataset["Randomised"])[1, ] == "no")
  
})

#Test inverse function#
test_that("Inverse function works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                     type = "variable", stat = "max", func = "inv", cmissing=FALSE)
  
  furthest = 2
  closest = 2
  duration  <- (furthest - closest) + 1
  maxmodno  <- (duration * (duration + 1))/2
  
  expect_true(is.list(test))
  expect_false(is.na((test[[1]]$BestModel)[1]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 2)
  
  expect_equal(length(which(is.na(test[[1]]$Dataset[, 4]))), 0)
  expect_true(test[[1]]$Dataset[, 8] == "inv")
  expect_true(ncol(test[[1]]$Dataset) >= 15)
  expect_equal(maxmodno, nrow(test[[1]]$Dataset))
  expect_true((test[[1]]$Dataset["Randomised"])[1, ] == "no")
  
})

##########################################################

#Test different cinterval values#

#Test cinterval = W
test_that("Weekly interval works", {
  
data(Mass, envir = environment())
data(MassClimate, envir = environment())

test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                   baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                   type = "variable", stat = "max", func = "lin",
                   cmissing=FALSE, cinterval = "week")

expect_true(is.list(test))
expect_false(is.na((test[[1]]$BestModel)[1]))

expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
expect_true(ncol(test[[1]]$BestModelData) >= 2)

})

#Test cinterval = M
test_that("Monthly interval works", {
  
data(Mass, envir = environment())
data(MassClimate, envir = environment())
  
test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                   baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                   type = "variable", stat = "max", func = "lin",
                   cmissing=FALSE, cinterval = "month")

expect_true(is.list(test))
expect_false(is.na((test[[1]]$BestModel)[1]))

expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
expect_true(ncol(test[[1]]$BestModelData) >= 2)
  
})

################################################################

#Error when you have NAs in the biological data
test_that("climatewin gives error when NAs are present in biological data", {
  
  data(MassClimate, envir = environment())
  Mass <- data.frame(Date = c("01/01/2014", "01/02/2014"), Mass = c(NA, 1))
  
  expect_error(climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                          baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                          type = "variable", stat = "max", func = "lin",
                          cmissing = FALSE, cinterval = "day"))  

  })

################################################################

#Test that cross validation works
test_that("Does cross validation work?", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- climatewin(xvar = list(MassClimate$Temp), cdate = MassClimate$Date, bdate = Mass$Date, 
                     baseline = lm(Mass ~ 1, data = Mass), furthest = 2, closest = 2, 
                     type = "variable", stat = "max", func = "lin",
                     cmissing = FALSE, cinterval = "day", cvk = 2)
  
  expect_true(is.list(test))
  expect_false(is.na((test[[1]]$BestModel)[1]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 2)
  
  expect_equal(test[[1]]$Dataset$K, 2)
  
})

################################################################

#Test mean centring
test_that("Mean centring is functioning", {
  
  data(Offspring, envir = environment())
  data(OffspringClimate, envir = environment())
  Offspring$Year <- lubridate::year(as.Date(Offspring$Date, format = "%d/%m/%Y"))
  
  test <- climatewin(xvar = list(OffspringClimate$Temp), cdate = OffspringClimate$Date, 
                     bdate = Offspring$Date, baseline = lm(Offspring ~ 1, data = Offspring), furthest = 2, closest = 2, 
                     type = "variable", stat = "max", func = "lin",
                     cmissing = FALSE, cinterval = "day", centre = Offspring$Year)
  
  expect_true(is.list(test))
  expect_false(is.na((test[[1]]$BestModel)[1]))
  
  expect_equal(length(which(is.na(test[[1]]$BestModelData))), 0)
  expect_true(ncol(test[[1]]$BestModelData) >= 3)
  
  expect_equal(length(which(is.na(test[[1]]$Dataset[, 4]))), 0)
  expect_equal(length(which(is.na(test[[1]]$Dataset[, 5]))), 0)
  expect_true(ncol(test[[1]]$Dataset) >= 14)
  expect_true((test[[1]]$Dataset["Randomised"])[1, ] == "no")
  
})