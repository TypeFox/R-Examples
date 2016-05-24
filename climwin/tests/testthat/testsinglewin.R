# Test singlewin function #

# Test that singlewin outputs SingleBestModel and SingleBestModelData
# Check that coefficients of SingleBestModel are not NA
# Check that there are no NAs in SingleBestModelData
# Check that SingleBestModelData has at least 2 columns

test_that("singlewin creates an output", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                    cdate = MassClimate$Date, bdate = Mass$Date,
                    baseline = lm(Mass$Mass~1), furthest = 72, closest = 15,
                    stat = "mean", func = "lin",
                    type = "variable", cmissing = FALSE, cinterval = "day")
  
  expect_true(is.list(test))  
  expect_false(is.na(test[[1]][1]))
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

test_that("cinterval W works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                    cdate = MassClimate$Date, bdate = Mass$Date,
                    baseline = lm(Mass ~ 1, data = Mass), furthest = 1, closest = 0,
                    stat = "mean", func = "lin",
                    type = "variable", cmissing = FALSE, cinterval = "week")
  
  expect_true(is.list(test))  
  expect_false(is.na(test[[1]][1]))
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

test_that("cinterval M works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                    cdate = MassClimate$Date, bdate = Mass$Date,
                    baseline = lm(Mass ~ 1, data = Mass), furthest = 1, closest = 0,
                    stat = "mean", func = "lin",
                    type = "variable", cmissing = FALSE, cinterval = "month")
  
  expect_true(is.list(test))  
  expect_false(is.na(test[[1]][1]))
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

###############################################################################################

# Test different settings of cmissing #

#When cmissing is TRUE and no NA is present#
test_that("No errors return when cmissing TRUE and full dataset", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                    cdate = MassClimate$Date, bdate = Mass$Date, 
                    baseline = lm(Mass ~ 1, data = Mass), 
                    furthest = 2, closest = 2, 
                    type = "variable", stat = "max", 
                    func = "lin", cmissing = TRUE)
  
  expect_true(is.list(test))  
  expect_false(is.na(test[[1]][1]))
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

#When cmissing is TRUE and NA is present#
test_that("No errors return when cmissing TRUE with NAs", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  MassClimate2 <- MassClimate[-491, ]
  test <- singlewin(xvar = list(Temp = MassClimate2$Temp), 
                    cdate = MassClimate2$Date, bdate = Mass$Date, 
                    baseline = lm(Mass ~ 1, data = Mass), 
                    furthest = 2, closest = 0, 
                    type = "variable", stat = "max", 
                    func = "lin", cmissing = TRUE)
  
  expect_true(is.list(test))  
  expect_false(is.na(test[[1]][1]))
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

#When cmissing is FALSE and NA is present#
test_that("No errors return when cmissing FALSE with NAs", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  MassClimate2 <- MassClimate[-491, ]
  expect_error(singlewin(xvar = list(Temp = MassClimate2$Temp), 
                         cdate = MassClimate2$Date, bdate = Mass$Date, 
                         baseline = lm(Mass ~ 1, data = Mass), 
                         furthest = 2, closest = 2, 
                         type = "variable", stat = "max", func = "lin", 
                         cmissing = FALSE))
  
})


##########################################################

# Test different types of models #

# Test glm models #
test_that("glm models can run", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                    cdate = MassClimate$Date, bdate = Mass$Date, 
                    baseline = glm(Mass ~ 1, data = Mass, family = poisson), 
                    furthest = 2, closest = 2, 
                    type = "variable", stat = "max", 
                    func = "lin", cmissing = FALSE)
  
  expect_true(is.list(test))
  expect_false(is.na((test[[1]])[1]))
  
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

test_that("lmer models can run", {
  
  data(Offspring, envir = environment())
  data(OffspringClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = OffspringClimate$Temp), 
                    cdate = OffspringClimate$Date, 
                    bdate = Offspring$Date, 
                    baseline = lmer(Offspring ~ 1 + (1|BirdID), data = Offspring),  
                    furthest = 2, closest = 2, type = "variable", 
                    stat = "max", func = "lin", cmissing = FALSE)
  
  expect_true(is.list(test))
  expect_false(is.na(fixef(test[[1]])[1]))
  expect_false(is.na(fixef(test[[1]])[2]))
  
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

##########################################################

# Test fixed and variable #

# Test fixed window#
test_that("Fixed window works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                    cdate = MassClimate$Date, bdate = Mass$Date, 
                    baseline = lm(Mass ~ 1, data = Mass), 
                    furthest = 2, closest = 2, 
                    type = "fixed", cutoff.day = 20, cutoff.month = 5, 
                    stat = "max", func = "lin", cmissing = FALSE)
  
  expect_true(is.list(test))
  
})

##########################################################

# Test slope stat #
test_that("slope stats work", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                    cdate = MassClimate$Date, bdate = Mass$Date, 
                    baseline = lm(Mass ~ 1, data = Mass), 
                    furthest = 2, closest = 1, 
                    type = "variable", stat = "slope", 
                    func = "lin", cmissing = FALSE)
  
  expect_true(is.list(test))  
  expect_false(is.na(test[[1]][1]))
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

##########################################################

#Test different functions for fitting climate#

#Test quadratic function#
test_that("Quadratic function works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                    cdate = MassClimate$Date, bdate = Mass$Date, 
                    baseline = lm(Mass ~ 1, data = Mass), 
                    furthest = 2, closest = 2, 
                    type = "variable", stat = "max", 
                    func = "quad", cmissing = FALSE)
  
  expect_true(is.list(test))  
  expect_false(is.na(test[[1]][1]))
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

#Test cubic function#
test_that("Cubic function works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                    cdate = MassClimate$Date, bdate = Mass$Date, 
                    baseline = lm(Mass ~ 1, data = Mass), 
                    furthest = 2, closest = 2, 
                    type = "variable", stat = "max", 
                    func = "cub", cmissing = FALSE)
  
  expect_true(is.list(test))  
  expect_false(is.na(test[[1]][1]))
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

#Test log function#
test_that("Log function works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                    cdate = MassClimate$Date, bdate = Mass$Date, 
                    baseline = lm(Mass ~ 1, data = Mass), 
                    furthest = 2, closest = 2, 
                    type = "variable", stat = "max", 
                    func = "log", cmissing = FALSE)
  
  expect_true(is.list(test))  
  expect_false(is.na(test[[1]][1]))
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

#Test log function#
test_that("Inverse function works", {
  
  data(Mass, envir = environment())
  data(MassClimate, envir = environment())
  
  test <- singlewin(xvar = list(Temp = MassClimate$Temp), 
                    cdate = MassClimate$Date, bdate = Mass$Date, 
                    baseline = lm(Mass ~ 1, data = Mass), 
                    furthest = 2, closest = 2, 
                    type = "variable", stat = "max", 
                    func = "inv", cmissing = FALSE)
  
  expect_true(is.list(test))  
  expect_false(is.na(test[[1]][1]))
  expect_equal(length(which(is.na(test[[2]]))), 0)
  expect_true(ncol(test[[2]]) >= 2)
  
})

################################################################

#Error when you have NAs in the biological data
test_that("singlewin gives error when NAs are present in biological data", {
  
  data(MassClimate, envir = environment())
  Mass <- data.frame(Date = c("01/01/2014", "01/02/2014"), Mass = c(NA, 1))
  
  expect_error(singlewin(xvar = list(Temp = MassClimate$Temp), 
                         cdate = MassClimate$Date, bdate = Mass$Date, 
                         baseline = lm(Mass ~ 1, data = Mass), 
                         furthest = 2, closest = 2, 
                         type = "variable", stat = "max", 
                         func = "lin", cmissing = FALSE))
  
})
