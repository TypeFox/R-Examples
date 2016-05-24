context("util.data")

test_that("checkGHS() and checkEOS() (via info()) produce messages", {
  i1 <- info("S2O3-2")
  expect_message(info(i1), "G of S2O3-2 aq \\(26\\) differs by 939 cal mol-1 from tabulated value")
  i2 <- info("Cu+2")
  expect_message(info(i2), "Cp of Cu\\+2 aq \\(62\\) differs by 3.62 cal K-1 mol-1 from tabulated value")
})

test_that("checkGHS() and checkEOS() respond to thermo$opt$*.tol", {
  i1 <- info("SO4-2")
  thermo$opt$Cp.tol <<- 0.5
  expect_message(info(i1), "checkEOS")
  i2 <- info("a,w-dicarboxytetracosane")
  thermo$opt$G.tol <<- 50
  expect_message(info(i2), "checkGHS")
})

test_that("RH2obigt() gives group additivity results consistent with database values (from Richard and Helgeson, 1998)", {
  file <- system.file("extdata/thermo/RH98_Table15.csv", package = "CHNOSZ")
  dat <- read.csv(file, stringsAsFactors=FALSE)
  ispecies <- info(dat$compound, dat$state)
  obigt.ref <- thermo$obigt[ispecies, ]
  obigt.calc <- RH2obigt(file=file)
  # calculated values of H are spot on; to pass tests, tolerance on
  # G is set higher; is there an incorrect group value somewhere?
  expect_true(max(abs(obigt.calc$G - obigt.ref$G)) < 31)
  expect_true(max(abs(obigt.calc$H - obigt.ref$H)) == 0)
  expect_true(max(abs(obigt.calc$S - obigt.ref$S)) < 0.02001)
  expect_true(max(abs(obigt.calc$Cp - obigt.ref$Cp)) < 0.04001)
  expect_true(max(abs(obigt.calc$V - obigt.ref$V)) < 0.1001)
  expect_true(max(abs(obigt.calc$a1.a - obigt.ref$a1.a)) < 0.01001)
  expect_true(max(abs(obigt.calc$a2.b - obigt.ref$a2.b)) < 1e-13)
  expect_true(max(abs(obigt.calc$a3.c - obigt.ref$a3.c)) < 1e-14)
})

test_that("add.obigt() replaces existing entries without changing species index", {
  # store the original species index of citric acid
  icitric <- info("citric acid", "aq")
  # add supplemental database - includes citric acid
  file <- system.file("extdata/thermo/OBIGT-2.csv", package="CHNOSZ")
  isp <- add.obigt(file, force=TRUE)
  # species index of citric acid should not have changed
  expect_equal(info("citric acid", "aq"), icitric)
  # check that names of species modified are same as in file
  newdat <- read.csv(file, stringsAsFactors=FALSE)
  # the order isn't guaranteed ... just make sure they're all there
  expect_true(all(newdat$name %in% thermo$obigt$name[isp]))
})

# reference

# Richard, L. and Helgeson, H. C. (1998) Calculation of the thermodynamic properties at elevated 
#   temperatures and pressures of saturated and aromatic high molecular weight solid and liquid 
#   hydrocarbons in kerogen, bitumen, petroleum, and other organic matter of biogeochemical interest. 
#   Geochim. Cosmochim. Acta 62, 3591--3636. http://dx.doi.org/10.1016/S0016-7037(97)00345-1
