# copyright (C) 2014-2016 A.Rebecq
library(testthat)

context("Test calibration functions on big dataset")

test_that("Calibration functions check out with Calmar", {
  
  data("population_test")
  
  sample <- dataPop[dataPop$weight > 0,]

  wCalLin <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                         , method="linear", description=FALSE)
  
  wCalRaking <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                         , method="raking", description=FALSE, pct=FALSE)
  
  wCalLogit <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                            , method="logit", bounds=c(0.2,1.3), description=FALSE)
  
  
  expect_equal(wCalLin, poptest_calmar$weight_cal_lin, tolerance=1e-6)
  expect_equal(wCalRaking, poptest_calmar$weight_cal_raking, tolerance=1e-6)
  expect_equal(wCalLogit, poptest_calmar$weight_cal_logit, tolerance=1e-6)
  
  popTotal <- 50000
  
  wCalLin2 <- calibration(data=sample, marginMatrix=table_margins_2, colWeights="weight"
                         , method="linear", description=FALSE, popTotal=popTotal, pct=TRUE)
  
  wCalRaking2 <- calibration(data=sample, marginMatrix=table_margins_2, colWeights="weight"
                            , method="raking", description=FALSE, popTotal=popTotal, pct=TRUE)
  
  wCalLogit2 <- calibration(data=sample, marginMatrix=table_margins_2, colWeights="weight"
                           , method="logit", bounds=c(0.2,1.3), description=FALSE, popTotal=popTotal, pct=TRUE)
  
  
  expect_equal(wCalLin2, poptest_calmar$weight_cal_lin_2, tolerance=1e-6)
  expect_equal(wCalRaking2, poptest_calmar$weight_cal_raking_2, tolerance=1e-6)
  expect_equal(wCalLogit2, poptest_calmar$weight_cal_logit_2, tolerance=1e-6)
  
  ## Test that calibrated estimators are equal to population totals
  totalsCalVar <- sapply(table_margins_1[,1], function(x) { return(sum(dataPop[,x])) })
  expect_equal( sapply(table_margins_1[,1], function(x) { return(HTtotal(sample[,x], wCalLin)) }),
                totalsCalVar,
                tolerance=1e-6)
  expect_equal( sapply(table_margins_1[,1], function(x) { return(HTtotal(sample[,x], wCalLogit)) }),
                totalsCalVar,
                tolerance=1e-6)
  expect_equal( sapply(table_margins_1[,1], function(x) { return(HTtotal(sample[,x], wCalRaking2)) }),
                totalsCalVar,
                tolerance=1e-6)
  
  ########### Tests with non-response and scale factor
  ## TODO : check parameters for logit calibration
  sampleNR <- sample[sample$responding==1,]
  
  wCalLinNR <- calibration(data=sampleNR, marginMatrix=table_margins_1, colWeights="weight"
                         , method="linear", description=FALSE, scale=TRUE)
  wCalRakingNR <- calibration(data=sampleNR, marginMatrix=table_margins_1, colWeights="weight"
                           , method="raking", description=FALSE, scale=TRUE)
  wCalLogitNR <- calibration(data=sampleNR, marginMatrix=table_margins_1, colWeights="weight"
                           , method="logit", bounds=c(0.1,2.0), description=FALSE, scale=TRUE, popTot=50000)
  
  expect_equal(wCalLinNR, poptest_calmar_nr$weight_cal_lin, tolerance=1e-6)
  expect_equal(wCalRakingNR, poptest_calmar_nr$weight_cal_raking, tolerance=1e-6)
  expect_equal(wCalLogitNR, poptest_calmar_nr$weight_cal_logit, tolerance=1e-6)
  
#   testDistrib <- poptest_calmar_nr$weight_cal_logit/sampleNR$weight * sum(sampleNR$weight) / 50000 
#   print(summary(testDistrib))
  
  wCalLinNR2 <- calibration(data=sampleNR, marginMatrix=table_margins_2, popTot = 50000, pct=T, colWeights="weight"
                           , method="linear", description=FALSE, scale=TRUE)
  wCalRakingNR2 <- calibration(data=sampleNR, marginMatrix=table_margins_2, popTot = 50000, pct=T, colWeights="weight"
                              , method="raking", description=FALSE, scale=TRUE)
  wCalLogitNR2 <- calibration(data=sampleNR, marginMatrix=table_margins_2, popTot = 50000, pct=T, colWeights="weight"
                             , method="logit", bounds=c(0.1,2.0), description=FALSE, scale=TRUE)
  
  expect_equal(wCalLinNR2, poptest_calmar_nr$weight_cal_lin_2, tolerance=1e-6)
  expect_equal(wCalRakingNR2, poptest_calmar_nr$weight_cal_raking_2, tolerance=1e-6)
  expect_equal(wCalLogitNR2, poptest_calmar_nr$weight_cal_logit_2, tolerance=1e-6)
  
})

test_that("Test margin stats", {
  
  data("population_test")
  
  sample <- dataPop[dataPop$weight > 0,]
  
  testStats1 <- marginStats(sample, table_margins_1, colWeights = "weight")
  testStats2 <- marginStats(sample, table_margins_2, colWeights = "weight", pct=T, popTotal = 50000)
  
  expect_equal(testStats1[13,1], 18.94)
  expect_equal(testStats1[14,2], 20.02)
  expect_equal(testStats1[15,3], 1.44)
  
  expect_equal(testStats2[17,3], 0.27)
  expect_equal(testStats2[18,1], 30.21)
  expect_equal(testStats2[16,2], 10)
  
  ## TODO : test marginStats with calibration weights
  ## (on penalized calibration for instance)
  sample$wCal <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                      , method="linear", description=FALSE, popTotal = 50000)
  
  testCosts <- rep(Inf, length(table_margins_1[,1]))

  sample$wCal_penal <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
                             , method="linear", description=FALSE, costs=testCosts, popTotal = 50000)
  
  testStats3 <- marginStats(sample, table_margins_1, colWeights = "weight", colCalibratedWeights = "wCal", popTotal = 50000)
  expect_equal(testStats3[,4], rep(0, length(testStats3[,1])))
  
  testStats4 <- marginStats(sample, table_margins_1, colWeights = "weight", colCalibratedWeights = "wCal_penal", popTotal = 50000)
  expect_equal(testStats4[,4], rep(0, length(testStats4[,1])))
  
  #### TODO : improve testing for penalized calib
#   testCosts2 <- rep(Inf, length(table_margins_1[,1]))
#   testCosts2[1] <- 1
#   testCosts2[2] <- 1
#   testCosts2[3] <- 100
#   testCosts2[4] <- 1
#   testCosts2[5] <- 1
#   testCosts2[2] <- 100
#   testCosts2[3] <- 1
#   
#   print(testCosts2)
  
  # TODO : set gap
#   sample$wCal_penal2 <- calibration(data=sample, marginMatrix=table_margins_1[c(1:3),], colWeights="weight"
#                                    , method="linear", description=TRUE, costs=testCosts2[c(1:3)], popTotal=50000, gap=0.1, uCostPenalized = 1e-3)
# 
#   sample$wCal_penal2 <- calibration(data=sample, marginMatrix=table_margins_1, colWeights="weight"
#                                     , method="linear", description=TRUE, costs=testCosts2, popTotal=50000, gap=0.2, uCostPenalized = 1e-4)
#   
  
  # testStats5 <- marginStats(sample, table_margins_1[c(1:3),], colWeights = "weight", colCalibratedWeights = "wCal_penal2", popTotal = 50000)
  # testStats5 <- marginStats(sample, table_margins_1, colWeights = "weight", colCalibratedWeights = "wCal_penal2", popTotal = 50000)

  # print(testStats5)
  
})