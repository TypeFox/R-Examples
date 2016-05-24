context("equilibrate")

# set up some simple systems
# one acid
basis("CHNOS+")
species("acetic acid")
aone <- suppressMessages(affinity())
# acids
species(c("formic acid", "formate", "acetate"))
aacid <- suppressMessages(affinity())
# acids plus a sulfur species
species("H2S")
aacidS <- suppressMessages(affinity())
# proteins
species(delete=TRUE)
species(c("LYSC_CHICK", "MYG_PHYCA", "RNAS1_BOVIN", "CYC_BOVIN"))
aprot <- suppressMessages(affinity())

test_that("equilibrate() gives expected messages and errors for balance calculation", {
  # the following error is triggered by equil.react, not equil.boltzmann
  expect_error(equilibrate(aone), "at least two species needed")
  expect_message(equilibrate(aacid), "balance: from moles of CO2")
  expect_message(equilibrate(aacid), "n.balance is 2 1 1 2 ")
  expect_message(equilibrate(aacid), "loga.balance is -2.221848")
  expect_message(equilibrate(aacid, loga.balance=-3), "loga.balance is -3")
  expect_error(equilibrate(aacid, balance="length"), "some species are not proteins")
  expect_message(equilibrate(aacidS), "balance: from numeric argument value")
  expect_message(equilibrate(aacidS), "n.balance is 1 1 1 1 1")
  expect_message(equilibrate(aacidS), "loga.balance is -2.301029")
  expect_error(equilibrate(aacidS, balance="CO2"), "some species have no CO2 in the formation reaction")
  expect_message(equilibrate(aprot), "balance: from protein length")
  expect_message(equilibrate(aprot), "n.balance is 129 153 124 104")
  expect_message(equilibrate(aprot), "loga.balance is -0.292429")
  expect_message(equilibrate(aprot, normalize=TRUE), "using 'normalize' for molar formulas")
})

test_that("equilibrate() gives expected messages and errors for species selection", {
  # an error if we select no species
  expect_error(equilibrate(aacid, ispecies=numeric()), "the length of ispecies is zero")
  # an error if all affinities are NA
  aNA <- aacid
  aNA$values[1:2] <- NA
  expect_error(equilibrate(aNA, ispecies=1:2), "all species have NA affinities")
  # a message if we select only certain of the species
  expect_message(equilibrate(aacid, ispecies=1:2), "using 2 of 4 species")
})

test_that("equilibrate() keeps the same total loga.balance for normalize=TRUE or FALSE", {
  # use the proteins
  e.norm <- equilibrate(aprot, normalize=TRUE)
  e <- equilibrate(aprot)
  # the total activity of the balance in the two cases
  sumact.balance.norm <- sum(10^unlist(e.norm$loga.equil)*e.norm$m.balance)
  sumact.balance <- sum(10^unlist(e$loga.equil)*e$n.balance)
  expect_equal(sumact.balance.norm, sumact.balance)
})

test_that("equilibrate() reproduces an example from the literature", {
  # the reference values are the equilibrium logarithms of activities
  # of sulfur species at logfO2=-30 from Seewald, 2001
  # we name them here because S5O6-2 isn't on the plot at logfO2=-30, 
  # and to get them in order
  species.ref <- c("S3O6-2", "S2O6-2", "S2O4-2", "S3-2", "S2-2", "S2O3-2", "HSO3-", "SO2", "HSO4-", "H2S")
  # these values were read from the plot using g3data (http://www.frantz.fi/software/g3data.php)
  loga.ref <- c(-28.82, -24.70, -22.10, -14.19, -12.12, -11.86, -8.40, -7.40, -6.54, -1.95)
  # set up the system - see ?diagram for an example showing the entire plot
  basis("CHNOS+")
  basis(c("pH", "O2"), c(5, -30))
  # we include here all the species shown by Seewald, 2001
  species(c("H2S", "S2-2", "S3-2", "S2O3-2", "S2O4-2", "S3O6-2", "S5O6-2", "S2O6-2", "HSO3-", "SO2", "HSO4-"))
  a <- affinity(T=325, P=350)
  # loga.balance=-2 signifies 10 mmolal total sulfur
  e <- equilibrate(a, loga.balance=-2)
  # get the calculated activities of the reference species
  loga.equil <- unlist(e$loga.equil[match(species.ref, e$species$name)])
  # the test... the tolerance may seem high, but consider that the reference values
  # were read from a plot with 30 logfO2 units spanning 4 inches
  expect_true(all(abs(loga.equil-loga.ref) < 0.36))
})

test_that("equilibrate() can be used for huge values of Astar", {
  ## working out some bugs and testing new 'method' argument 20151109

  ## first, demonstrate that equil.reaction works where equil.boltzmann doesn't
  # minimal example: Astar=c(0, 0), n.balance=c(1, 1), loga.balance=0
  # results in equal activities of two species
  eb0 <- equil.boltzmann(c(0, 0), c(1, 1), 0)
  expect_equal(unlist(eb0), rep(log10(0.5), 2))
  # Astar=c(-330, -330)
  # result is NaN (we probably get an Inf-Inf somewhere)
  eb330 <- equil.boltzmann(c(-330, -330), c(1, 1), 0)
  expect_equal(unlist(eb330), rep(NaN, 2))
  # (fixed bug: while loop in equil.reaction tested a NaN value)
  # (dlogadiff.dAbar <- 0 / 0)
  er330 <- equil.reaction(c(-330, -330), c(1, 1), 0)
  expect_equal(er330, eb0)

  ## second, set up extreme test case and show boltzmann method produces NaN (is.na)
  basis("CHNOS")
  basis("O2", 200)
  species(c("glycine", "alanine", "proline"))
  a <- affinity()
  expect_message(eb <- equilibrate(a, balance=1), "using boltzmann method")
  expect_true(all(is.na(unlist(eb$loga.equil))))

  ## third, check we can use method="reaction"
  expect_message(er1 <- equilibrate(a, balance=1, method="reaction"), "using reaction method")
  expect_false(any(is.na(unlist(er1$loga.equil))))
  # is it an equilibrium solution?
  species(1:3, unlist(er1$loga.equil))
  a1 <- affinity()
  expect_equal(diff(range(unlist(a1$values))), 0)

  ## third, check that we can use arbitrary numeric balance specification
  # (balance <> 1 here means equilibrate will call equil.reaction)
  expect_message(er11 <- equilibrate(a, balance=1.000001), "using reaction method")
  species(1:3, unlist(er11$loga.equil))
  a11 <- affinity()
  expect_equal(unlist(a1$values), unlist(a11$values))

  ## fourth, check that equil.boltzmann won't run for balance <> 1
  expect_error(equilibrate(a, balance=1.000001, method="boltzmann"), "won't run equil.boltzmann")
})


# references

# Seewald, J. S. (2001) 
#   Aqueous geochemistry of low molecular weight hydrocarbons at elevated temperatures and
#   pressures: Constraints from mineral buffered laboratory experiments
#   Geochim. Cosmochim. Acta 65, 1641--1664. http://dx.doi.org/10.1016/S0016-7037(01)00544-0
