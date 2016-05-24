context("subcrt")

# delete the basis definition in case there is one
basis(delete=TRUE)

test_that("unbalanced reactions give a warning", {
  expect_warning(subcrt(c("glucose", "ethanol"), c(-1, 3)), "reaction was unbalanced, missing H-6O3")
})

test_that("unbalanced reactions are balanced given sufficient basis species", {
  basis("CHNOS")
  # since it doesn't alter the species indices of the basis species, this can come second ...
  add.obigt()
  s <- subcrt(c("malic acid", "citric acid"), c(-1, 1))
  expect_equal(s$reaction$coeff, c(-1, 1, -2, -1, 1.5))
  expect_equal(s$reaction$name, c("malic acid", "citric acid", "CO2", "water", "oxygen"))
})

test_that("phase transitions of minerals give expected messages and results", {
  iacanthite <- info("acanthite", "cr2")
  expect_message(subcrt(iacanthite), "subcrt: some points below transition temperature for acanthite cr2 \\(using NA for G\\)")
  expect_message(subcrt(iacanthite), "subcrt: some points above transition temperature for acanthite cr2 \\(using NA for G\\)")
  expect_equal(subcrt("acanthite")$out$acanthite$state, c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3))
})

test_that("heat capacity of minerals are consistent with literature", {
  # from Helgeson et al., 1978 Table 6
  # here, set P to NA so that water() is not called to get Psat
  expect_equal(subcrt("wollastonite", T=c(200, 500, 800), P=NA)$out$wollastonite$Cp, c(25.4, 28.3, 29.9), tolerance=1e-2)
})

test_that("standard Gibbs energies of reactions involving aqueous species are consistent with the literature", {
  # choose units of Joules for subcrt()
  E.units("J")
  # from Amend and Shock, 2001 Table 3
  T <- c(2, 18, 25, 37, 45, 55, 70, 85, 100, 115, 150, 200)
  # standard Gibbs energies in kJ/mol
  DG0.H2O <- c(78.25, 79.34, 79.89, 80.90, 81.63, 82.59, 84.13, 85.78, 87.55, 89.42, 94.22, 102.21)
  # H2O(l) = H+ + OH-
  sout.H2O <- subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T=T)$out
  # from Amend and Shock, 2001 Table 4.3 Reaction A1
  DG0.A1 <- c(-263.94, -263.45, -263.17, -262.62, -262.20, -261.63, -260.67, -259.60, -258.44, -257.18, -253.90, -248.44)
  # H2(aq) + 0.5O2(aq) = H2O(l)
  sout.A1 <- subcrt(c("H2", "O2", "H2O"), "aq", c(-1, -0.5, 1), T=T)$out
  # from Amend and Shock, 2001 Table 5.1 NO(g) and NO(aq)
  DG0.NO.g <- c(91.39, 88.04, 86.57, 84.03, 82.33, 80.20, 76.99, 73.75, 70.50, 67.23, 59.53, 48.38)
  DG0.NO.aq <- c(104.56, 102.87, 102.06, 100.58, 99.53, 98.15, 95.98, 93.67, 91.25, 88.72, 82.46, 72.77)
  # NO(g) = NO(aq)  ... greater tolerance, our values for NO(aq) differ slightly from AS01
  sout.NO <- subcrt("NO", c("gas", "aq"), c(-1, 1), T=T)$out
  # from Amend and Shock, 2001 Table 5.1 Reaction B10
  DG0.B10 <- c(-268.85, -268.01, -267.50, -266.46, -265.66, -264.55, -262.66, -260.54, -258.19, -255.63, -248.81, -237.06)
  # NH3(aq) + 1.5O2(aq) = H+ + NO2- + H2O(l)
  sout.B10 <- subcrt(c("NH3", "O2", "H+", "NO2-", "H2O"), "aq", c(-1, -1.5, 1, 1, 1), T=T)$out
  # from Amend and Shock, 2001 Table 6.3 Reaction C7
  DG0.C7 <- c(-1695.30, -1686.90, -1682.80, -1675.30, -1670.00, -1663.10, -1652.00, -1640.30, -1628.00, -1615.20, -1583.50, -1533.00)
  # 5S2O3-2 + H2O(l) + 4O2(aq) = 6SO4-2 + 2H+ + 4S(s)
  s.C7 <- subcrt(c("S2O3-2", "H2O", "O2", "SO4-2", "H+", "S"), c("aq", "liq", "aq", "aq", "aq", "cr"), c(-5, -1, -4, 6, 2, 4), T=T)
  sout.C7 <- s.C7$out
  # from Amend and Shock, 2001 Table 8.3 Reaction E12
  DG0.E12 <- c(132.52, 132.26, 132.29, 132.49, 132.74, 133.15, 133.98, 135.04, 136.31, 137.79, 141.97, 149.53)
  # 4(2-)propanol(aq) + 3CO2(aq) + 2H2O(l) = 3CH4(aq) + 4lactic acid(aq)
  sout.E12 <- subcrt(c("2-propanol", "CO2", "H2O", "CH4", "lactic acid"), c(-4, -3, -2, 3, 4), T=T)$out
  # now the tests, tolerances set to lowest order of magnitute to pass
  expect_equal(sout.H2O$G/1000, DG0.H2O, tolerance=1e-4)
  expect_equal(sout.A1$G/1000, DG0.A1, tolerance=1e-4)
  # greater tolerance, our values for NO(aq) differ slightly from AS01
  expect_equal(sout.NO$G/1000, DG0.NO.aq - DG0.NO.g, tolerance=1e-2)
  expect_equal(sout.B10$G/1000, DG0.B10, tolerance=1e-3)
  # we can check that sulfur has expected phase transitions
  expect_equal(s.C7$state$sulfur, c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3))
  expect_equal(sout.C7$G/1000, DG0.C7, tolerance=1e-4)
  # this one is on hold until the HKF parameters of 2-propanol can be located
  #expect_equal(sout.E12$G/1000, DG0.E12, 1e-4)
  # also todo: COS(g) in Amend and Helgeson, 2001 Table 7.2?
  # return to our favourite units
  E.units("cal")
})

test_that("subzero degree C calculations are possible", {
  ## start with H2O
  s.H2O <- subcrt("H2O", T=c(-20.1, seq(-20, 0)), P=1)$out$water
  # we shouldn't get anything at -20.1 deg C
  expect_equal(s.H2O$G[1], NA_real_)
  # we should get something at -20 deg C
  expect_equal(floor(s.H2O$G[2]), -56001)
  # for historical reasons, an input temperature of 0 was converted to 0.01
  expect_equal(s.H2O$T[22], 0.01)
})

test_that("calculations using IAPWS-95 are possible", {
  thermo$opt$water <<- "IAPWS95"
  sb <- subcrt(c("H2O", "Na+"), T=c(-30, -20, 0, 10), P=1)$out
  # the test is not a demanding numerical comparison, more that we got numbers and no error
  expect_that(all(sb$`Na+`$G < sb$water$G), is_true())
  # clean up
  suppressMessages(data(thermo))
})

# references

# Amend, J. P. and Shock, E. L. (2001) 
#   Energetics of overall metabolic reactions of thermophilic and hyperthermophilic Archaea and Bacteria. 
#   FEMS Microbiol. Rev. 25, 175--243. http://dx.doi.org/10.1016/S0168-6445(00)00062-0

# Helgeson, H. C., Delany, J. M., Nesbitt, H. W. and Bird, D. K. (1978) 
#   Summary and critique of the thermodynamic properties of rock-forming minerals. 
#   Am. J. Sci. 278-A, 1--229. http://www.worldcat.org/oclc/13594862

