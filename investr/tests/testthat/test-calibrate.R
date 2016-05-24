# The following tests are for calibration with simple linear regression models.
context("Inverse estimation in the SLR model")


test_that("output matches answers to Graybill and Iyer (1996, chap. 6)", {
    
  # Thermostat example from Graybill and Iyer (1996, p. 431)
  thermom <- data.frame(temp = seq(from = 96, to = 110, by = 2), 
                        read = c(95.71, 98.16, 99.52, 102.09, 103.79, 106.18, 
                                 108.14, 110.21))
  thermom.cal1 <- calibrate(thermom, y0 = 104)
  thermom.cal2 <- calibrate(thermom, y0 = 100, level = 0.9)
  expect_that(round(thermom.cal1$estimate, 3), equals(103.995))
  expect_that(round(thermom.cal1$lower, 1), equals(103.4))
  expect_that(round(thermom.cal1$upper, 1), equals(104.6))
  expect_that(round(thermom.cal2$estimate, 1), equals(100.1))
  expect_that(round(thermom.cal2$lower, 2), equals(99.63))
  expect_that(round(thermom.cal2$upper, 2), equals(100.59))
  
  # Reaction chamber example from Graybill and Iyer (1996, p. 433)
  chamber <- data.frame(dial = seq(from = 0, to = 100, by = 10), 
                        temp = c(206.36, 225.52, 252.18, 289.33, 318.11, 349.49, 
                                 383.03, 410.70, 444.40, 469.14, 501.16))
  chamber.reg <- calibrate(chamber, y0 = 400, mean.response = TRUE, 
                           level = 0.99)
  
  # Expectations for regulation
  expect_that(round(chamber.reg$estimate, 1), equals(66.5))
  expect_that(round(chamber.reg$lower, 2), equals(65.07))
  expect_that(round(chamber.reg$upper, 2), equals(68.03))
  
  # Crystal weight example from Graybill and Iyer (1996, p. 434)
  crystal.lm <- lm(weight ~ time, data = crystal)
  crystal.reg <- calibrate(crystal.lm, y0 = 5, mean.response = TRUE, 
                           level = 0.9)
  
  # Expectations for calibration
  expect_equal(round(crystal.reg$estimate, 2), 9.93)
  expect_equal(round(crystal.reg$lower, 2), 8.65)
  expect_equal(round(crystal.reg$upper, 2), 11.05)
  
})


test_that("errors are handled appropriately", {
  
  # Simulated data
  set.seed(101)
  x <- rep(seq(from = 0, to = 10, length = 10), 2)
  y <-  3 + 0.01*x + rnorm(length(x), sd = 0.5)
  d1 <- data.frame(x, y)
  d2 <- list(x = 1:11, y = 1:10 + rnorm(10, sd = 1))
  
  # Expectations
  expect_that(calibrate(d1, y0 = 3), gives_warning())
  expect_that(calibrate(d1, y0 = 3, mean.response = TRUE), gives_warning())
  expect_that(calibrate(d1, y0 = 2), throws_error())
  expect_that(calibrate(d1, y0 = 2.5, mean.response = TRUE), throws_error())
  expect_that(calibrate(d2, y0 = 2.5), throws_error())
  expect_that(calibrate(y ~ x + I(x^2), y0 = 2.5), throws_error())
  
})


test_that("approximate standard error is correct", {
  
  # Crystal weight example from Graybill and Iyer (1996, p. 434)
  crystal.lm <- lm(weight ~ time, data = crystal)
  crystal.cal <- calibrate(crystal.lm, y0 = 5, interval = "Wald")
  crystal.reg <- calibrate(crystal.lm, y0 = 5, interval = "Wald", 
                           mean.response = TRUE)
  
  # Calculate and compare standard error using invest and car::deltaMethod
#   covmat.cal <- diag(3)
#   covmat.cal[1:2, 1:2] <- vcov(crystal.lm)
#   covmat.cal[3, 3] <- summary(crystal.lm)$sigma^2
#   coefs <- unname(coef(crystal.lm))
#   params <- c(b0 = coefs[1], b1 = coefs[2], y0 = 5)
  se.cal <- 2.211698 #car::deltaMethod(params, g = "(y0-b0)/b1", vcov. = covmat.cal)$SE
  se.reg <- 0.6658998 #car::deltaMethod(crystal.lm, g = "(5-b0)/b1", 
                             #parameterNames = c("b0", "b1"))$SE

  # Expectations
  expect_that(crystal.cal$se, equals(se.cal, tol = 1e-04)) # small diff
  expect_that(crystal.reg$se, equals(se.reg, tol = 1e-04))
  
})


test_that("all methods produce equivalent results", {
  
  # Generate some data
  set.seed(101)
  x <- rep(1:10, each = 3)
  y <- 2 + 3 * x + rnorm(length(x), sd = 1)
  d <- data.frame(x = x, y = y)
  
  # Matrix method - inversion interval
  cal1 <- calibrate(cbind(x, y), y0 = 15, mean.response = FALSE)
  
  # data.frame method - inversion interval
  cal2 <- calibrate(data.frame(x, y), y0 = 15, mean.response = FALSE)
  
  # list method - inversion interval
  cal3 <- calibrate(list(x, y), y0 = 15, mean.response = FALSE)
  
  # formula method - inversion interval
  cal4 <- calibrate(y ~ x, y0 = 15, mean.response = FALSE)
  
  # formula method w/ data - inversion interval
  cal5 <- calibrate(y ~ x, data = d, y0 = 15, mean.response = FALSE)
  
  # lm method - inversion interval
  cal6 <- calibrate(lm(y ~ x), y0 = 15, mean.response = FALSE)
  
  # formula method w/ transformations - inversion interval
  cal7 <- calibrate(exp(log(y)) ~ sqrt(x^2), y0 = 15, mean.response = FALSE)

  # Matrix method - Wald interval
  cal8 <- calibrate(cbind(x, y), y0 = 15, mean.response = FALSE, 
                    interval = "Wald")
  
  # Matrix method - Wald interval
  cal9 <- calibrate(data.frame(x, y), y0 = 15, mean.response = FALSE, 
                    interval = "Wald")
  
  # data.frame method - Wald interval
  cal10 <- calibrate(list(x, y), y0 = 15, mean.response = FALSE, 
                     interval = "Wald")
  
  # formula method - Wald interval
  cal11 <- calibrate(y ~ x, y0 = 15, mean.response = FALSE, interval = "Wald")
  
  # formula method w/ data - Wald interval
  cal12 <- calibrate(y ~ x, data = d, y0 = 15, mean.response = FALSE, 
                    interval = "Wald")
  
  # lm method - Wald interval
  cal13 <- calibrate(lm(y ~ x), y0 = 15, mean.response = FALSE, 
                     interval = "Wald")
  
  # formula method method w/transformations - Wald interval
  cal14 <- calibrate(exp(log(y)) ~ sqrt(x^2), y0 = 15, mean.response = FALSE, 
                     interval = "Wald")
  
  # These should all be identical
  expect_identical(cal1, cal2)
  expect_identical(cal1, cal3)
  expect_identical(cal1, cal4)
  expect_identical(cal1, cal5)
  expect_identical(cal1, cal6)
  expect_equal(cal1, cal7)  # Why are these two not identical?
  
  # These should all be identical
  expect_identical(cal8, cal9)
  expect_identical(cal8, cal10)
  expect_identical(cal8, cal11)
  expect_identical(cal8, cal12)
  expect_identical(cal8, cal13)
  expect_equal(cal8, cal14)  # Why are these two not identical?
  
})


test_that("errors get handled apprropriately", {
  
  # Nonlinear least squares fit
  nls.fit <- nls(weight ~ theta1/(1 + exp(theta2 + theta3 * log(conc))),
                 start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
                 data = nasturtium)
  
  # Multiple linear regression
  mlr.fit1 <- lm(weight ~ time + I(time ^ 2), data = crystal)
  mlr.fit2 <- lm(cbind(weight, weight ^ 2) ~ time, data = crystal)
  
  # Expectations
  expect_error(calibrate(nls.fit, y0 = c(309, 296, 419)))
  expect_error(calibrate(mlr.fit1, y0 = c(309, 296, 419)))
  expect_error(calibrate(mlr.fit1, y0 = c(309, 296, 419)))
  
})


test_that("multiple inference procedures work", {
  
  # Crystal weight example from Graybill and Iyer (1996, p. 434)
  crystal.lm <- lm(weight ~ time, data = crystal)
  crystal.cal <- calibrate(crystal.lm, y0 = 5, interval = "Wald")
  crystal.cal.multi <- calibrate(crystal.lm, y0 = 5, interval = "Wald", 
                                 adjust = "Scheffe", k = 1)
  
  # Expectations
  expect_equal(crystal.cal, crystal.cal.multi)
  
})
