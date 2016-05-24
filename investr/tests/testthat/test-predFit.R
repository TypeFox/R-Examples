# The following tests are for testing the prediction function and estimated
# standard errors of the fitted values
context("Prediction")

# DNase data from the dataframes package
DNase1 <- data.frame(conc = c(0.04882812, 0.04882812, 0.19531250, 0.19531250, 
                              0.39062500, 0.39062500, 0.78125000, 0.78125000, 
                              1.56250000, 1.56250000, 3.12500000, 3.12500000, 
                              6.25000000, 6.25000000, 12.50000000, 12.50000000),
                     density = c(0.017, 0.018, 0.121, 0.124, 0.206, 0.215, 
                                 0.377, 0.374, 0.614, 0.609, 1.019, 1.001, 
                                 1.334, 1.364, 1.730, 1.710))

# Nonlinear model fit
DNase1_nls <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
                  data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
                    
test_that("predFit matches results from PROC NLIN in SAS", {
  
  # New data
  DNase1.new <- data.frame(conc = c(8.0, 12.0, 1.0, 10.0, 5.5))
  
  # Predictions, standard errors, etc.
  DNase1.conf <- predFit(DNase1_nls, se.fit = TRUE, interval = "confidence")
  DNase1.conf2 <- predFit(DNase1_nls, se.fit = TRUE, newdata = DNase1.new, 
                          interval = "confidence")
  DNase1.pred <- predFit(DNase1_nls, se.fit = TRUE, interval = "prediction")
  DNase1.pred2 <- predFit(DNase1_nls, se.fit = TRUE, newdata = DNase1.new, 
                          interval = "prediction")
  
  # Fitted value standard errors from PROC NLIN in SAS/STATS
  SAS.se.fit <- c(0.002899, 0.002899, 0.006079, 0.006079, 0.007461, 0.007461, 
                  0.007704, 0.007704, 0.007581, 0.007581, 0.009277, 0.009277,
                  0.009038, 0.009038, 0.012925, 0.012925)
  SAS.se.fit.new <- c(0.008857, 0.012272, 0.007548, 0.010014, 0.009314)
  
  # Confidence limits from PROC NLIN in SAS/STATS
  SAS.conf.lwr <- c(0.02442, 0.02442, 0.09892, 0.09892, 0.19246, 0.19246,
                    0.35768, 0.35768, 0.61640, 0.61640, 0.96082, 0.96082,
                    1.34799, 1.34799, 1.68706, 1.68706)
  SAS.conf.upr <- c(0.03694, 0.03694, 0.12518, 0.12518, 0.22470, 0.22470,
                    0.39097, 0.39097, 0.64916, 0.64916, 1.00090, 1.00090,
                    1.38704, 1.38704, 1.74291, 1.74291)
  SAS.conf.lwr.new <- c(1.48029, 1.67025, 0.43872, 1.58988, 1.27678)
  SAS.conf.upr.new <- c(1.51856, 1.72327, 0.47133, 1.63315, 1.31703)
  
  # Prediction limits from PROC NLIN in SAS/STATS
  SAS.pred.lwr <- c(-0.01126, -0.01126,  0.06855,  0.06855,  0.16409,  0.16409,  
                    0.32965,  0.32965,  0.58819,  0.58819,  0.93481,  0.93481,
                    1.32168,  1.32168,  1.66500,  1.66500)
  SAS.pred.upr <- c(0.07262, 0.07262, 0.15555, 0.15555, 0.25307, 0.25307,
                    0.41901, 0.41901, 0.67736, 0.67736, 1.02692, 1.02692,
                    1.41335, 1.41335, 1.76498, 1.76498)
  SAS.pred.lwr.new <- c(1.45376, 1.64754, 0.41047, 1.56474, 1.25081)
  SAS.pred.upr.new <- c(1.54510, 1.74598, 0.49958, 1.65829, 1.34300)
  
  # Expectations for original data
  expect_true(all.equal(round(DNase1.conf$se.fit, 6), SAS.se.fit))
  expect_true(all.equal(round(DNase1.pred$se.fit, 6), SAS.se.fit))
  expect_true(all.equal(round(DNase1.conf$fit[, "lwr"], 5), SAS.conf.lwr, tol = 1e-05))
  expect_true(all.equal(round(DNase1.conf$fit[, "upr"], 5), SAS.conf.upr, tol = 1e-05))
  expect_true(all.equal(round(DNase1.pred$fit[, "lwr"], 5), SAS.pred.lwr, tol = 1e-05))
  expect_true(all.equal(round(DNase1.pred$fit[, "upr"], 5), SAS.pred.upr, tol = 1e-05))
  
  # Expectations for new data
  expect_true(all.equal(round(DNase1.conf2$se.fit, 6), SAS.se.fit.new))
  expect_true(all.equal(round(DNase1.pred2$se.fit, 6), SAS.se.fit.new))
  expect_true(all.equal(round(DNase1.conf2$fit[, "lwr"], 5), SAS.conf.lwr.new, tol = 1e-05))
  expect_true(all.equal(round(DNase1.conf2$fit[, "upr"], 5), SAS.conf.upr.new, tol = 1e-05))
  expect_true(all.equal(round(DNase1.pred2$fit[, "lwr"], 5), SAS.pred.lwr.new, tol = 1e-05))
  expect_true(all.equal(round(DNase1.pred2$fit[, "upr"], 5), SAS.pred.upr.new, tol = 1e-05))
  
})


test_that("predFit works properly on 'special' nls fits", {
  
  # Using conditional linearity
  DNase1_nls_2 <- nls(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
                      data = DNase1,
                      start = list(xmid = 0, scal = 1),
                      algorithm = "plinear")
  
  # Using selfStart
  DNase1_nls_3 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), 
                      data = DNase1)
  
  
  # Using Port's nl2sol algorithm
  DNase1_nls_4 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
                      data = DNase1,
                      start = list(Asym = 3, xmid = 0, scal = 1),
                      algorithm = "port")
  
  # Predictions
  pred <- predFit(DNase1_nls, se.fit = TRUE, interval = "prediction")
  pred_3 <- predFit(DNase1_nls_3, se.fit = TRUE, interval = "prediction")
  pred_4 <- predFit(DNase1_nls_4, se.fit = TRUE, interval = "prediction")
  
  # Expectations
  expect_error(predFit(DNase1_nls_2))
  expect_true(all.equal(pred$fit[, "fit"], pred_3$fit[, "fit"], tol = 1e-06))
  expect_true(all.equal(pred$fit[, "fit"], pred_4$fit[, "fit"], tol = 1e-06))
  expect_true(all.equal(pred$fit[, "lwr"], pred_3$fit[, "lwr"], tol = 1e-06))
  expect_true(all.equal(pred$fit[, "lwr"], pred_4$fit[, "lwr"], tol = 1e-06))
  expect_true(all.equal(pred$fit[, "upr"], pred_3$fit[, "upr"], tol = 1e-06))
  expect_true(all.equal(pred$fit[, "upr"], pred_4$fit[, "upr"], tol = 1e-06))
  expect_true(all.equal(pred$se.fit, pred_3$se.fit, tol = 1e-06))
  expect_true(all.equal(pred$se.fit, pred_4$se.fit, tol = 1e-06))
  
})


test_that("predFit matches output from stats::predict", {

  # Simulate some data
  set.seed(101)  # for reproducibilty
  x <- rep(1:10, each = 3)
  y <- 1 + 2 * x + rnorm(length(x), sd = 1)
  d <- data.frame("x" = x, "y" = y)

  # Fit some linear models
  lm1 <- lm(y ~ x, data = d)
  lm2 <- lm(y ~ x)

  # Predictions only
  pred.investr <- predFit(lm1)
  pred.stats <- predict(lm1)

  # Predictions and confidence intervals
  pred.investr.conf <- predFit(lm1, interval = "confidence")
  pred.stats.conf <- predict(lm1, interval = "confidence")

  # Predictions and prediction intervals
  pred.investr.pred <- predFit(lm1, interval = "prediction")
  pred.stats.pred <- predict(lm1, interval = "prediction")

  # Predictions and standard errors
  pred.investr.se <- predFit(lm1, se.fit = TRUE)
  pred.stats.se <- predict(lm1, se.fit = TRUE)

  # Predictions, confidence intervals, and standard errors
  pred.investr.se.conf <- predFit(lm1, se.fit = TRUE, interval = "confidence")
  pred.stats.se.conf <- predict(lm1, se.fit = TRUE, interval = "confidence")

  # Predictions, prediction intervals, and standard errors
  pred.investr.se.pred <- predFit(lm1, se.fit = TRUE, interval = "prediction")
  pred.stats.se.pred <- predict(lm1, se.fit = TRUE, interval = "prediction")

  # Expectations
  expect_equal(pred.investr, pred.stats)
  expect_equal(pred.investr.conf, pred.stats.conf)
  expect_equal(pred.investr.pred, pred.stats.pred)
  expect_equal(pred.investr.se$se.fit, pred.stats.se$se.fit)
  expect_equal(pred.investr.se.conf$se.fit, pred.stats.se$se.fit)
  expect_equal(pred.investr.se.pred$se.fit, pred.stats.se$se.fit)
  expect_equal(pred.investr.se.conf$fit, pred.stats.se.conf$fit)
  expect_equal(pred.investr.se.pred$fit, pred.stats.se.pred$fit)
  expect_equal(predFit(lm1, se.fit = TRUE, interval = "prediction"),
               predFit(lm2, se.fit = TRUE, interval = "prediction"))
  
})


context("Working-Hotelling band")


test_that("predFit reproduces example from pplied Linear Statistical Models", {
  
  # Toluca company data from Applied Linear Statistical Models (ALSM)
  lotsize <- c(80, 30, 50, 90, 70, 60, 120, 80, 100, 50, 40, 70, 90, 
               20, 110, 100, 30, 50, 90, 110, 30, 90, 40, 80, 70)
  workhrs <- c(399, 121, 221, 376, 361, 224, 546, 352, 353, 157, 160, 252, 389, 
               113, 435, 420, 212, 268, 377, 421, 273, 468, 244, 342, 323)
  toluca <- data.frame(lotsize, workhrs)
  
  # Simple linear regression model
  toluca_lm <- lm(workhrs ~ lotsize, data = toluca)
  toluca_nls <- nls(workhrs ~ b0 + b1*lotsize, data = toluca,
                    start = list(b0 = 62.37, b1 = 3.57))
  
  
  # Working-Hotelling band
  wh_lm <- predFit(toluca_lm, newdata = data.frame("lotsize" = c(30, 65, 100)), 
                   se.fit = TRUE, interval = "confidence", level = 0.9,
                   adjust = "Scheffe")
  wh_nls <- predFit(toluca_nls, newdata = data.frame("lotsize" = c(30, 65, 100)), 
                    se.fit = TRUE, interval = "confidence", level = 0.9,
                    adjust = "Scheffe")
  # Answer from ALSM
  ans <- cbind("fit" = c(169.5, 294.4, 419.4),
               "lwr" = c(131.2, 272.0, 387.2),
               "upr" = c(207.8, 316.8,451.6))
  rownames(ans) <- 1:3
  
  # Expectations
  expect_equal(wh_lm, wh_nls, tol = 1e-05, check.attributes = FALSE)
  expect_identical(unname(round(wh_lm$se.fit, digits = 2)), c(16.97, 9.92, 14.27))
  expect_identical(round(wh_lm$fit, digits = 1), ans)
  
})