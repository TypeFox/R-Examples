context("affinity")

# clear out any previous basis definition or database alterations
suppressMessages(data(thermo))

test_that("errors come as expected, and output gives T and P in user's units", {
  expect_error(affinity(iprotein=7), "basis species are not defined")
  basis("CHNOS")
  expect_error(affinity(), "species have not been defined")
  species("5a(H),14b(H)-cholestane")
  a.C_bar <- affinity(T=c(0, 100, 10), P=c(10, 1000, 10))
  expect_equal(range(a.C_bar$vals[[1]]), c(0, 100))
  expect_equal(range(a.C_bar$vals[[2]]), c(10, 1000))
  T.units("K")
  P.units("MPa")
  a.K_MPa <- affinity(T=c(273.15, 373.15, 10), P=c(1, 100, 10))
  expect_equal(range(a.K_MPa$vals[[1]]), c(273.15, 373.15))
  expect_equal(range(a.K_MPa$vals[[2]]), c(1, 100))
  # different units, same T,P ... same affinities
  expect_equal(a.C_bar$values, a.K_MPa$values)
  # go back to original units for the remaining tests
  T.units("C")
  P.units("bar")
})

test_that("pe, pH and Eh are correctly handled", {
  basis("CHNOSe")
  species(c("HS-", "H2S", "SO4-2"))
  Eh <- c(-1, 1)
  pe <- convert(Eh, "pe", T=convert(100, "K"))
  a.Eh <- affinity(Eh=Eh, T=100)
  a.pe <- affinity(pe=pe, T=100)
  # they should give the same result
  # ... except for names(dim(.)), so set check.attributes=FALSE
  expect_equal(a.Eh$values, a.pe$values, check.attributes=FALSE)
  # the variables should have the right names
  expect_equal(c(a.Eh$vars, a.pe$vars), c("Eh", "pe"))
  # now for an Eh-pH example
  pH <- c(0, 14)
  a <- affinity(pH=pH, Eh=Eh)
  expect_equal(a$vars, c("pH", "Eh"))
  expect_equal(range(a$vals[[1]]), pH)
  expect_equal(range(a$vals[[2]]), Eh)
  expect_equal(length(a$vals[[2]]), 128)
  # since Eh has to be reconstructed, check it's done correctly
  a129 <- affinity(pH=pH, Eh=c(Eh, 129))
  expect_equal(length(a129$vals[[2]]), 129)

  ## a transect of hotter, more oxidizing and more acidic
  ## has not been working since at least 0.9-7
  ##T <- c(25, 50, 100, 125, 150)
  ##Eh <- c(-1, -0.5, 0, 0.5, 1)
  ##pH <- c(10, 8, 6, 4, 2)
  ##a <- affinity(T=T, Eh=Eh, pH=pH)
})

test_that("affinity() in 3D returns values consistent with manual calculation", {
  # our "manual" calculation will be for H2(aq) + 0.5O2(aq) = H2O(l)
  # the equilibrium constants at 25 and 100 degrees C
  # (the logK are tested against literature values in test-subcrt.R)
  logK.25 <- subcrt(c("H2", "O2", "H2O"), "aq", c(-1, -0.5, 1), T=25)$out$logK
  logK.100 <- subcrt(c("H2", "O2", "H2O"), "aq", c(-1, -0.5, 1), T=100)$out$logK
  # the value of A/2.303RT at 25 degrees and logaH2=-10, logaO2=-10 and logaH2O=0
  A.2303RT.25.10.10 <- logK.25 - ( (-1)*(-10) + (-0.5)*(-10) )
  # the value of A/2.303RT at 100 degrees and logaH2=-5, logaO2=-10 and logaH2O=0
  A.2303RT.100.5.10 <- logK.100 - ( (-1)*(-5) + (-0.5)*(-10) )
  # set up basis and species
  basis(c("H2", "O2"), "aq")
  species("H2O")
  # we will run affinity() in 3D
  # T = 0, 25, 50, 75, 100, 125 degrees
  # log_a(H2) = -20, -15, -10, -5, 0
  # log_a(O2) = -20, -15, -10, -5, 0
  # first test: the dimensions are correct
  a.logK <- affinity(T=c(0, 125, 6), H2=c(-20, 0, 5), O2=c(-20, 0, 5), property="logK")
  expect_equal(dim(a.logK$values[[1]]), c(6, 5, 5), check.names=FALSE)
  # second and third tests: the logK values used by affinity() are correct
  expect_equal(a.logK$values[[1]][2, 3, 3], logK.25)
  expect_equal(a.logK$values[[1]][5, 4, 3], logK.100)
  # fourth and fifth tests: the A/2.303RT values returned by affinity() are correct
  a.A <- affinity(T=c(0, 125, 6), H2=c(-20, 0, 5), O2=c(-20, 0, 5))
  expect_equal(a.A$values[[1]][2, 3, 3], A.2303RT.25.10.10)
  expect_equal(a.A$values[[1]][5, 4, 3], A.2303RT.100.5.10)
})

test_that("'iprotein' gives consistent results on a transect", {
  # from Dick and Shock, 2011, values of A/2.303RT for the per-residue
  # formation reactions of overall model proteins at five sampling sites
  # at Bison Pool, with different temperature, pH and log_a(H2)
  # these are the maximum values for each site from Table 5 in the paper
  A.2303RT_ref <- c(-18.720, -27.894, -35.276, -36.657, -41.888)
  # the measured temperatures and pHs
  T <- c(93.3, 79.4, 67.5, 65.3, 57.1)
  pH <- c(7.350, 7.678, 7.933, 7.995, 8.257)
  # Eq. 24 of the paper
  H2 <- -11+T*3/40
  # remove "RESIDUE" entries in thermo$obigt (clutter from first test)
  data(thermo)
  basis(c("HCO3-", "H2O", "NH3", "HS-", "H2", "H+"),
    "aq", c(-3, 0, -4, -7, 999, 999))
  sites <- c("N", "S", "R", "Q", "P")
  aa <- read.aa(system.file("extdata/protein/DS11.csv", package="CHNOSZ"))
  ip <- add.protein(aa[1:5, ])
  # to reproduce, we need use the "old" parameters for [Met]
  add.obigt()
  a <- affinity(T=T, pH=pH, H2=H2, iprotein=ip)
  # divide A/2.303RT by protein length
  pl <- protein.length(ip)
  A.2303RT <- t(sapply(a$values, c)) / pl
  # find the maximum for each site
  A.2303RT_max <- apply(A.2303RT, 2, max)
  # we're off a bit in the second decimal ... 
  # maybe becuase of rounding of the aa composition?
  expect_equal(A.2303RT_max, A.2303RT_ref, tolerance=1e-3)
  # todo: add comparison with results from loading proteins via species()
})

test_that("affinity() for proteins (with/without 'iprotein') returns same value as in previous package versions", {
  # our test case is CSG_HALJP because it has no methionine
  # (aqueous [Met] was updated in 0.9-9)
  # these values were calculated using versions 0.6, 0.8 and 0.9-7
  # (25 degrees C, 1 bar, basis species "CHNOS" or "CHNOS+")
  A.2303RT.nonionized <- -3795.297
  A.2303RT.ionized <- -3075.222
  # first for nonionized protein
  basis("CHNOS")
  # try it with iprotein
  ip <- iprotein("CSG_HALJP")
  expect_equal(affinity(iprotein=ip)$values[[1]][1], A.2303RT.nonionized, tolerance=1e-6)
  # then with the protein loaded as a species
  species("CSG_HALJP")
  expect_equal(affinity()$values[[1]][1], A.2303RT.nonionized, tolerance=1e-6)
  # now for ionized protein
  basis("CHNOS+")
  expect_equal(affinity(iprotein=ip)$values[[1]][1], A.2303RT.ionized, tolerance=1e-6)
  species("CSG_HALJP")
  expect_equal(affinity()$values[[1]][1], A.2303RT.ionized, tolerance=1e-6)
})

test_that("affinity() for proteins keeps track of pH on 2-D calculations", {
  # (relates to the "thisperm" construction in A.ionization() )
  basis("CHNOS+")
  species("LYSC_CHICK")
  a1 <- affinity(pH=c(6, 8, 3))
  a2 <- affinity(pH=c(6, 8, 3), T=c(0, 75, 4))
  expect_equal(as.numeric(a1$values[[1]]), a2$values[[1]][, 2])
})

test_that("IS can be constant or variable", {
  # inspired by an error from demo("phosphate")
  # > a25 <- affinity(IS=c(0, 0.14), T=T[1])
  # ...
  # Error in subcrt(species = c(1017L, 20L, 19L), property = "logK", T = 298.15,  : 
  #   formal argument "IS" matched by multiple actual arguments
  basis("CHNOPS+")
  species(c("PO4-3", "HPO4-2", "H2PO4-"))
  a0 <- affinity()
  a1 <- affinity(IS=0.14)
  a2 <- affinity(IS=c(0, 0.14))
  expect_equal(unlist(lapply(a2$values, head, 1)), unlist(a0$values))
  expect_equal(unlist(lapply(a2$values, tail, 1)), unlist(a1$values))
})
