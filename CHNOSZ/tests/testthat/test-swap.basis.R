context("swap.basis")

# clear out any previous basis definition or database alterations
suppressMessages(data(thermo))

test_that("swap.basis raises errors when needed", {
  expect_error(swap.basis(c("CO2", "H2O")), "requires an existing basis definition")
  basis("CHNOS+")
  expect_error(swap.basis(c("CO2", "H2O")), "two species must be identified")
  expect_error(swap.basis(c("CO2", "H2O"), c("HCO3-", "H2O")), "can only swap one species for one species")
  expect_error(swap.basis("CH4", "C2H5OH"), "basis species .* is not defined")
  expect_error(swap.basis("CO2", "C60"), "is not available")
  expect_error(swap.basis("CO2", "H2"), "overdetermined")
})

test_that("basis.logact only accepts defined elements", {
  # setup basis species with two elements: C and H
  basis(c("graphite", "H2"), c("cr", "gas"))
  # we can't get basis activities with one element
  expect_error(basis.logact(c(C=1)), "number of elements in 'emu' is less than those in basis")
  # get some potentials of C, H and O
  ispecies <- info(c("ethane", "propane", "acetic acid", "propanoic acid"))  
  w <- run.wjd(ispecies, as.chemical.formula(colMeans(i2A(ispecies))))
  ep <- equil.potentials(w)
  # try to calculate log activities of basis species: get an error
  expect_error(basis.logact(ep), "element\\(s\\) O not found in basis")
})

test_that("equil.potentials - basis.logact - element.mu makes a roundtrip at 25 and 100 degrees C", {
  basis(c("graphite", "H2", "O2"), c("cr", "gas", "gas"))
  ispecies <- info(c("ethane", "propane", "acetic acid", "propanoic acid"))  
  # at 25 degrees C
  w25 <- run.wjd(ispecies, as.chemical.formula(colMeans(i2A(ispecies))))
  ep25 <- equil.potentials(w25)
  bl25 <- basis.logact(ep25)
  # set the activities of the bais species
  basis(names(bl25), bl25)
  # element.mu() calculates the chemical potentials of the elements from the current setting of basis species
  expect_equal(element.mu(), ep25)
  # at 100 degrees C
  w100 <- run.wjd(ispecies, as.chemical.formula(colMeans(i2A(ispecies))), T=100)
  ep100 <- equil.potentials(w100)
  bl100 <- basis.logact(ep100, T=100)
  basis(names(bl100), bl100)
  expect_equal(element.mu(T=100), ep100)
})
