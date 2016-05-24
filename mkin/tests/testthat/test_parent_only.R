# Copyright (C) 2015 Johannes Ranke
# Contact: jranke@uni-bremen.de

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>

context("Fitting of parent only models")

calc_dev.percent <- function(fitlist, reference, endpoints = TRUE, round_results = NULL) {
  dev.percent <- list()
  for (i in 1:length(fitlist)) {
    fit <- fitlist[[i]]
    if (endpoints) {
      results <- c(fit$bparms.optim, 
                   endpoints(fit)$distimes$DT50,
                   endpoints(fit)$distimes$DT90)
    } else {
      results <- fit$bparms.optim
    }
    if (!missing(round_results)) results <- round(results, round_results)
    dev.percent[[i]] <- abs(100 * ((reference - results)/reference))
  }
  return(dev.percent)
}

SFO <- mkinmod(parent = list(type = "SFO"))
FOMC <- mkinmod(parent = list(type = "FOMC"))
DFOP <- mkinmod(parent = list(type = "DFOP"))
HS <- mkinmod(parent = list(type = "HS"))
SFORB <- mkinmod(parent = list(type = "SFORB"))

test_that("Fits for FOCUS A deviate less than 0.1% from median of values from FOCUS report", {
  fit.A.SFO <- list(mkinfit("SFO", FOCUS_2006_A, quiet = TRUE))

  median.A.SFO <- as.numeric(lapply(subset(FOCUS_2006_SFO_ref_A_to_F, 
                                        dataset == "A", 
                                        c(M0, k, DT50, DT90)), "median"))

  dev.percent.A.SFO <- calc_dev.percent(fit.A.SFO, median.A.SFO)
  expect_equivalent(dev.percent.A.SFO[[1]] < 0.1, rep(TRUE, 4))

  # Fitting FOCUS A with FOMC is possible, but the correlation between
  # alpha and beta, when obtained, is 1.0000, and the fit sometimes failed on
  # Windows, as the Port algorithm did not converge (winbuilder, 2015-05-15)
  fit.A.FOMC <- try(list(mkinfit("FOMC", FOCUS_2006_A, quiet = TRUE)))
  if (!inherits(fit.A.FOMC, "try-error")) {

    median.A.FOMC <- as.numeric(lapply(subset(FOCUS_2006_FOMC_ref_A_to_F, 
                                          dataset == "A", 
                                          c(M0, alpha, beta, DT50, DT90)), "median"))

    dev.percent.A.FOMC <- calc_dev.percent(fit.A.FOMC, median.A.FOMC)
    # alpha and are beta ill-determined, do not compare those
    expect_equivalent(dev.percent.A.FOMC[[1]][c(1, 4, 5)] < 0.1, rep(TRUE, 3))
  }

  fit.A.DFOP <- list(mkinfit("DFOP", FOCUS_2006_A, quiet = TRUE))

  median.A.DFOP <- as.numeric(lapply(subset(FOCUS_2006_DFOP_ref_A_to_B, 
                                        dataset == "A", 
                                        c(M0, k1, k2, f, DT50, DT90)), "median"))

  dev.percent.A.DFOP <- calc_dev.percent(fit.A.DFOP, median.A.DFOP)
  #expect_equivalent(dev.percent.A.DFOP[[1]] < 0.1, rep(TRUE, 6)) # g/f is ill-determined
  expect_equivalent(dev.percent.A.DFOP[[1]][c(1, 2, 3, 5, 6)] < 0.1, rep(TRUE, 5))

  fit.A.HS <- list(mkinfit("HS", FOCUS_2006_A, quiet = TRUE))

  median.A.HS <- as.numeric(lapply(subset(FOCUS_2006_HS_ref_A_to_F, 
                                        dataset == "A", 
                                        c(M0, k1, k2, tb, DT50, DT90)), "median"))

  dev.percent.A.HS <- calc_dev.percent(fit.A.HS, median.A.HS)
  expect_equivalent(dev.percent.A.HS[[1]] < 0.1, rep(TRUE, 6))
})

test_that("Fits for FOCUS B deviate less than 0.1% from median of values from FOCUS report", {
  fit.B.SFO <- list(mkinfit("SFO", FOCUS_2006_B, quiet = TRUE))

  median.B.SFO <- as.numeric(lapply(subset(FOCUS_2006_SFO_ref_A_to_F, 
                                        dataset == "B", 
                                        c(M0, k, DT50, DT90)), "median"))

  dev.percent.B.SFO <- calc_dev.percent(fit.B.SFO, median.B.SFO)
  expect_equivalent(dev.percent.B.SFO[[1]] < 0.1, rep(TRUE, 4))

  fit.B.FOMC <- list(mkinfit("FOMC", FOCUS_2006_B, quiet = TRUE))

  median.B.FOMC <- as.numeric(lapply(subset(FOCUS_2006_FOMC_ref_A_to_F, 
                                        dataset == "B", 
                                        c(M0, alpha, beta, DT50, DT90)), "median"))

  dev.percent.B.FOMC <- calc_dev.percent(fit.B.FOMC, median.B.FOMC)
  expect_equivalent(dev.percent.B.FOMC[[1]] < 0.1, rep(TRUE, 5))

  fit.B.DFOP <- list(mkinfit("DFOP", FOCUS_2006_B, quiet = TRUE))

  median.B.DFOP <- as.numeric(lapply(subset(FOCUS_2006_DFOP_ref_A_to_B, 
                                        dataset == "B", 
                                        c(M0, k1, k2, f, DT50, DT90)), "median"))

  dev.percent.B.DFOP <- calc_dev.percent(fit.B.DFOP, median.B.DFOP)
  #expect_equivalent(dev.percent.B.DFOP[[1]] < 0.1, rep(TRUE, 6)) # g/f is ill-determined
  expect_equivalent(dev.percent.B.DFOP[[1]][c(1, 2, 3, 5, 6)] < 0.1, rep(TRUE, 5))

  fit.B.HS <- list(mkinfit("HS", FOCUS_2006_B, quiet = TRUE))

  median.B.HS <- as.numeric(lapply(subset(FOCUS_2006_HS_ref_A_to_F, 
                                        dataset == "B", 
                                        c(M0, k1, k2, tb, DT50, DT90)), 
                                   "median", na.rm = TRUE))

  dev.percent.B.HS <- calc_dev.percent(fit.B.HS, median.B.HS)
  expect_equivalent(dev.percent.B.HS[[1]] < 0.1, rep(TRUE, 6))

  fit.B.SFORB <- list(mkinfit(SFORB, FOCUS_2006_B, quiet=TRUE))
  dev.percent.B.SFORB <- calc_dev.percent(fit.B.SFORB, median.B.DFOP)
  expect_equivalent(dev.percent.B.SFORB[[1]][c(1, 5, 6)] < 0.1, rep(TRUE, 3))
})

test_that("Fits for FOCUS C deviate less than 0.1% from median of values from FOCUS report", {
  fit.C.SFO <- list(mkinfit("SFO", FOCUS_2006_C, quiet = TRUE))

  median.C.SFO <- as.numeric(lapply(subset(FOCUS_2006_SFO_ref_A_to_F, 
                                        dataset == "C", 
                                        c(M0, k, DT50, DT90)), "median"))

  dev.percent.C.SFO <- calc_dev.percent(fit.C.SFO, median.C.SFO)
  expect_equivalent(dev.percent.C.SFO[[1]] < 0.1, rep(TRUE, 4))

  fit.C.FOMC <- list(mkinfit("FOMC", FOCUS_2006_C, quiet = TRUE))

  median.C.FOMC <- as.numeric(lapply(subset(FOCUS_2006_FOMC_ref_A_to_F, 
                                        dataset == "C", 
                                        c(M0, alpha, beta, DT50, DT90)), "median"))

  dev.percent.C.FOMC <- calc_dev.percent(fit.C.FOMC, median.C.FOMC, 
                                         round_results = 2) # Not enough precision in FOCUS results
  expect_equivalent(dev.percent.C.FOMC[[1]] < 0.1, rep(TRUE, 5))

  fit.C.HS <- list(mkinfit("HS", FOCUS_2006_C, quiet = TRUE))

  median.C.HS <- as.numeric(lapply(subset(FOCUS_2006_HS_ref_A_to_F, 
                                        dataset == "C", 
                                        c(M0, k1, k2, tb, DT50, DT90)), "median"))

  dev.percent.C.HS <- calc_dev.percent(fit.C.HS, median.C.HS, round_results = c(2, 4, 6, 2, 2))
  # Not enouth precision in k2 available
  expect_equivalent(dev.percent.C.HS[[1]] < c(0.1, 0.1, 0.3, 0.1, 0.1, 0.1), rep(TRUE, 6))
})

test_that("SFO fits give approximately (0.001%) equal results with different solution methods", {
  fit.A.SFO.default <- mkinfit("SFO", FOCUS_2006_A, quiet = TRUE)$bparms.optim

  fits.A.SFO <- list()
  fits.A.SFO[[1]] <- mkinfit(SFO, FOCUS_2006_A, quiet = TRUE)
  fits.A.SFO[[2]] <- mkinfit(SFO, FOCUS_2006_A, quiet = TRUE, solution_type = "eigen")
  fits.A.SFO[[3]] <- mkinfit(SFO, FOCUS_2006_A, quiet = TRUE, solution_type = "deSolve")

  dev.percent <- calc_dev.percent(fits.A.SFO, fit.A.SFO.default, endpoints = FALSE)
  expect_equivalent(dev.percent[[1]] < 0.001, rep(TRUE, 2))
  expect_equivalent(dev.percent[[2]] < 0.001, rep(TRUE, 2))
  expect_equivalent(dev.percent[[3]] < 0.001, rep(TRUE, 2))
})

test_that("FOMC fits give approximately (0.001%) equal results with different solution methods", {
  fit.C.FOMC.default <- mkinfit("FOMC", FOCUS_2006_C, quiet = TRUE)$bparms.optim

  fits.C.FOMC <- list()
  fits.C.FOMC[[1]] <- mkinfit(FOMC, FOCUS_2006_C, quiet = TRUE)
  fits.C.FOMC[[2]] <- mkinfit(FOMC, FOCUS_2006_C, quiet = TRUE, solution_type = "deSolve")

  dev.percent <- calc_dev.percent(fits.C.FOMC, fit.C.FOMC.default, endpoints = FALSE)
  expect_equivalent(dev.percent[[1]] < 0.001, rep(TRUE, 3))
  expect_equivalent(dev.percent[[2]] < 0.001, rep(TRUE, 3))
})

test_that("DFOP fits give approximately (0.001%) equal results with different solution methods", {
  fit.C.DFOP.default <- mkinfit("DFOP", FOCUS_2006_C, quiet = TRUE)$bparms.optim

  fits.C.DFOP <- list()
  fits.C.DFOP[[1]] <- mkinfit(DFOP, FOCUS_2006_C, quiet = TRUE)
  fits.C.DFOP[[2]] <- mkinfit(DFOP, FOCUS_2006_C, quiet = TRUE, solution_type = "deSolve")

  dev.percent <- calc_dev.percent(fits.C.DFOP, fit.C.DFOP.default, endpoints = FALSE)
  expect_equivalent(dev.percent[[1]] < 0.001, rep(TRUE, 4))
  expect_equivalent(dev.percent[[2]] < 0.001, rep(TRUE, 4))
})

test_that("SFORB fits give approximately (0.002%) equal results with different solution methods", {
  fit.B.SFORB.default <- mkinfit(SFORB, FOCUS_2006_B, quiet=TRUE)$bparms.optim

  fits.B.SFORB <- list()          
  fits.B.SFORB[[1]] <- mkinfit(SFORB, FOCUS_2006_B, quiet=TRUE, solution_type = "eigen")
  fits.B.SFORB[[2]] <- mkinfit(SFORB, FOCUS_2006_B, quiet=TRUE, solution_type = "deSolve")
  dev.percent <- calc_dev.percent(fits.B.SFORB, fit.B.SFORB.default, endpoints = FALSE)
  expect_equivalent(dev.percent[[1]] < 0.001, rep(TRUE, 4))
  expect_equivalent(dev.percent[[2]] < 0.002, rep(TRUE, 4))
})
