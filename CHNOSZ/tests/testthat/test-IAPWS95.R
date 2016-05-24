context("IAPWS95")

test_that("calculations of Helmholtz free energy and its derivatives are consistent with reference cases", {
  ## reference values of these terms are listed Table 6.6 of Wagner and Pruss, 2002  
  p <- c("phi", "phi.delta", "phi.delta.delta", "phi.tau", "phi.tau.tau", "phi.delta.tau")
  # reference values for case 1: T=500 K, rho=838.025 kg m-3
  idealgas.ref.1 <- c(0.204797734e1, 0.384236747, -0.147637878, 0.904611106e1, -0.193249185e1, 0)
  residual.ref.1 <- c(-0.342693206e1, -0.364366650, 0.856063701, -0.581403435e1, -0.223440737e1, -0.112176915e1)
  # reference values for case 1: T=647 K, rho=358 kg m-3
  idealgas.ref.2 <- c(-0.156319605e1, 0.899441341, -0.808994726, 0.980343918e1, -0.343316334e1, 0)
  residual.ref.2 <- c(-0.121202657e1, -0.714012024, 0.475730696, -0.321722501e1, -0.996029507e1, -0.133214720e1)
  ## set up the problem
  # critical point constants
  T.critical <- 647.096 # K
  rho.critical <- 322 # kg m-3
  # T and rho for cases 1 and 2
  T <- c(500, 647)
  rho <- c(838.025, 358)
  # delta and tau for cases 1 and 2
  delta <- rho / rho.critical
  tau <- T.critical / T
  # calculated ideal gas and residual parts for case 1
  idealgas.calc.1 <- sapply(p, IAPWS95.idealgas, delta[1], tau[1])
  residual.calc.1 <- sapply(p, IAPWS95.residual, delta[1], tau[1])
  # calculated ideal gas and residual parts for case 2
  idealgas.calc.2 <- sapply(p, IAPWS95.idealgas, delta[2], tau[2])
  residual.calc.2 <- sapply(p, IAPWS95.residual, delta[2], tau[2])
  ## perform the tests
  # we almost get away without increasing the tolerance in any test ...
  expect_equal(idealgas.calc.1, idealgas.ref.1, check.attributes=FALSE)
  expect_equal(residual.calc.1, residual.ref.1, check.attributes=FALSE)
  expect_equal(idealgas.calc.2, idealgas.ref.2, check.attributes=FALSE)
  # ... however an offset is apparent in the value of the residual phi.delta.delta for case 2
  expect_equal(residual.calc.2, residual.ref.2, check.attributes=FALSE, tolerance=1e-5)
})

test_that("calculations of thermodynamic properties are consistent with reference values", {
  ## these are the properties we test - from Table 13.1 of Wagner and Pruss, 2002
  ## (speed of sound omitted as it's not currently implemented)
  p <- c("P", "H", "S", "Cv", "Cp")
  ## a selection of T and rho at vapor-liquid boundary
  ## (NOTE: excluding triple and critical points; we have some unresolved issues there)
  T <- c(274, 320, 368, 416, 464, 512, 560, 608, 647)
  rho.liquid <- c(999.843, 989.387, 961.984, 923.577, 875.125, 814.982, 737.831, 626.74, 357.34)
  rho.vapor <- c(0.00514, 0.07166, 0.50231, 2.1203, 6.5107, 16.409, 37.147, 84.173, 286.51)
  ## reference values
  P.ref <- c(0.000650, 0.010546, 0.084142, 0.39166, 1.2788, 3.2798, 7.1062, 13.681, 22.038)
  H.liquid.ref <- c(3.544, 196.170, 397.457, 601.396, 811.225, 1032.06, 1273.11, 1558.42, 2029.44)
  H.vapor.ref <- c(2502.46, 2585.71, 2667.37, 2737.09, 2785.91, 2803.05, 2771.24, 2646.01, 2148.56)
  S.liquid.ref <- c(0.0130, 0.6629, 1.2487, 1.7687, 2.2436, 2.6914, 3.1319, 3.6034, 4.3224)
  S.vapor.ref <- c(9.1331, 8.1302, 7.4169, 6.9025, 6.4994, 6.1504, 5.8071, 5.3922, 4.5065)
  cv.liquid.ref <- c(4.2156, 4.0416, 3.7950, 3.5560, 3.3523, 3.1884, 3.0722, 3.0617, 6.2344)
  cv.vapor.ref <- c(1.4191, 1.4627, 1.5428, 1.7138, 2.0015, 2.3697, 2.8225, 3.4538, 6.2740)
  cp.liquid.ref <- c(4.2171, 4.1807, 4.2101, 4.2892, 4.4513, 4.7616, 5.4239, 7.6169, 3905.2)
  cp.vapor.ref <- c(1.8852, 1.9417, 2.0601, 2.3334, 2.8561, 3.7263, 5.4099, 10.805, 5334.1)
  ## calculated values
  liquid.calc <- IAPWS95(p, T, rho.liquid)
  vapor.calc <- IAPWS95(p, T, rho.vapor)
  ## the tests
  # take P to 5 significant digits but not more than 6 decimals
  expect_equal(round(signif(liquid.calc$p, 5), 6), P.ref, tolerance=1e-3)
  expect_equal(round(signif(vapor.calc$p, 5), 6), P.ref, tolerance=1e-4)
  # take H to 6 significant digits but not more than 3 decimals
  expect_equal(round(signif(liquid.calc$h, 6), 3), H.liquid.ref, tolerance=1e-5)
  expect_that(round(signif(vapor.calc$h, 6), 3), equals(H.vapor.ref))  # spot on!
  # round S to 4 decimals
  expect_that(round(liquid.calc$s, 4), equals(S.liquid.ref))  # spot on!
  expect_equal(round(vapor.calc$s, 4), S.vapor.ref, tolerance=1e-4)
  # round cv to 4 decimals
  expect_equal(round(liquid.calc$cv, 4), cv.liquid.ref, tolerance=1e-4)
  expect_equal(round(vapor.calc$cv, 4), cv.vapor.ref, tolerance=1e-4)
  # take cp to 5 significant digits but not more than 4 decimals
  # note high tolerance setting: the highest temperature is the challenge
  expect_equal(round(signif(liquid.calc$cp, 5), 4), cp.liquid.ref, tolerance=1e0)
  expect_equal(round(signif(vapor.calc$cp, 5), 4), cp.vapor.ref, tolerance=1e-1)
})

test_that("calculations are possible at low temperatures", {
  # the sequences start at the lowest whole-number temperature (K) 
  # where the function returns a value of density at the given pressure
  expect_that(any(is.na(water.IAPWS95("rho", T=seq(234, 274, 3), P=rep(1, 14)))), is_false())
  expect_that(any(is.na(water.IAPWS95("rho", T=seq(227, 272, 3), P=rep(1000, 16)))), is_false())
})

# reference

# Wagner, W. and Pruss, A. (2002)
#   The IAPWS formulation 1995 for the thermodynamic properties of
#   ordinary water substance for general and scientific use
#   J. Phys. Chem. Ref. Data 31, 387--535. http://dx.doi.org/10.1063/1.1461829
