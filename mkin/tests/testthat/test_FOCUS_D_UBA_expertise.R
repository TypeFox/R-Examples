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

context("Results for FOCUS D established in expertise for UBA (Ranke 2014)")

SFO_SFO <- mkinmod(parent = list(type = "SFO", to = "m1"),
                   m1 = list(type = "SFO"), quiet = TRUE)
SFO_SFO.ff <- mkinmod(parent = list(type = "SFO", to = "m1"),
                      m1 = list(type = "SFO"), 
                      use_of_ff = "max", quiet = TRUE)

fit.default <- mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE)
fit.ff <- mkinfit(SFO_SFO.ff, FOCUS_2006_D, quiet = TRUE)

# Results are from p. 40

test_that("Fitted parameters are correct for FOCUS D", {
  expect_equivalent(round(fit.ff$bparms.optim, c(2, 4, 4, 4)), 
                    c(99.60, 0.0987, 0.0053, 0.5145))
})

test_that("FOCUS chi2 error levels are correct for FOCUS D", {
  expect_equivalent(round(100 * mkinerrmin(fit.ff)$err.min, 2), 
                    c(6.40, 6.46, 4.69))
})

test_that("DT50/90 are correct for FOCUS D when using formation fractions", {
  expect_equal(round(as.numeric(endpoints(fit.ff)$distimes["parent", ]), 2), 
               c(7.02, 23.33))
  expect_equal(round(as.numeric(endpoints(fit.ff)$distimes["m1", ]), 1), 
               c(131.8, 437.7))
})

test_that("DT50/90 are correct for FOCUS D when not using formation fractions", {
  expect_equal(round(as.numeric(endpoints(fit.default)$distimes["parent", ]), 2), 
               c(7.02, 23.33))
  expect_equal(round(as.numeric(endpoints(fit.default)$distimes["m1", ]), 1), 
               c(131.8, 437.7))
})

# References:
# Ranke (2014) Prüfung und Validierung von Modellierungssoftware als Alternative
# zu ModelMaker 4.0, Umweltbundesamt Projektnummer 27452

context("The t-test for significant difference from zero")

test_that("The t-value for fits using internal transformations corresponds with result from FME", {

  expect_equal(signif(summary(fit.default)$bpar[, "t value"], 5),
               c(parent_0 = 61.720, k_parent_sink = 12.777, k_parent_m1 = 24.248, k_m1_sink = 7.3486))

})

m_synth_DFOP_par.minff <- mkinmod(parent = mkinsub("DFOP", c("M1", "M2")),
                           M1 = mkinsub("SFO"),
                           M2 = mkinsub("SFO"),
                           use_of_ff = "min", quiet = TRUE)

fit_DFOP_par_c_2 <- mkinfit(m_synth_DFOP_par.minff, 
                          synthetic_data_for_UBA_2014[[12]]$data,
                          quiet = TRUE)

test_that("The t-value for fits using internal transformations corresponds with results from FME, synthetic data", {

  # Note that the k1 and k2 are exchanged in the untransformed fit evaluated with FME for this test
  expect_equal(signif(summary(fit_DFOP_par_c_2)$bpar[1:7, "t value"], 5),
               c(parent_0 = 80.054, k_M1_sink = 12.291, k_M2_sink = 10.588, 
                 f_parent_to_M1 = 21.4960, f_parent_to_M2 = 24.0890,
                 k1 = 16.1450, k2 = 8.1747))

})
