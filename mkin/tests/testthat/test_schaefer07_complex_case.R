# Copyright (C) 2014-2015 Johannes Ranke
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

# This test was migrated from a RUnit test inst/unitTests/runit.mkinfit.R

context("Complex test case from Schaefer et al. (2007) Piacenza paper")

schaefer07_complex_model <- mkinmod(
  parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
  A1 = list(type = "SFO", to = "A2"),
  B1 = list(type = "SFO"),
  C1 = list(type = "SFO"),
  A2 = list(type = "SFO"), use_of_ff = "max", quiet = TRUE)

schaefer07_long <- mkin_wide_to_long(schaefer07_complex_case, time = "time")

fit.default <- mkinfit(schaefer07_complex_model, schaefer07_long, quiet = TRUE)

test_that("Complex test case from Schaefer (2007) can be reproduced (10% tolerance)", {
  
  s <- summary(fit.default)
  r <- schaefer07_complex_results

  with(as.list(fit.default$bparms.optim), {
    r$mkin <<- c(
      k_parent,
      s$distimes["parent", "DT50"],
      s$ff["parent_A1"],
      k_A1,
      s$distimes["A1", "DT50"],
      s$ff["parent_B1"],
      k_B1,
      s$distimes["B1", "DT50"],
      s$ff["parent_C1"],
      k_C1,
      s$distimes["C1", "DT50"],
      s$ff["A1_A2"],
      k_A2,
      s$distimes["A2", "DT50"])
    }
  )
  r$means <- (r$KinGUI + r$ModelMaker)/2
  r$mkin.deviation <- abs(round(100 * ((r$mkin - r$means)/r$means), digits=1))
  expect_equal(r$mkin.deviation < 10, rep(TRUE, 14))
})

test_that("We avoid the local minumum with default settings", {
  # If we use optimisation algorithm 'Marq' we get a local minimum with a
  # sum of squared residuals of 273.3707
  # When using 'Marq', we need to give a good starting estimate e.g. for k_A2 in
  # order to get the optimum with sum of squared residuals 240.5686
  expect_equal(round(fit.default$ssr, 4), 240.5686)
})
