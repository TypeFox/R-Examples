context("basis")

# clear out any previous basis definition or database alterations
suppressMessages(data(thermo))

test_that("invalid basis definitions cause an error", {
  expect_error(basis(character()), "argument is empty")
  expect_error(basis(c("CO2", "CO2")), "names are not unique")
  expect_error(basis(c("CO2", "H2O")), "underdetermined")
  expect_error(basis(c("H2O", "O2", "H2")), "overdetermined")
  expect_error(basis(c("HCN", "H2O", "O2", "H2")), "singular")
  expect_error(basis(c("CN", "H2O", "O2", "H2")), "species not available")
  expect_error(basis(c("CN")), "species not available")
  expect_error(basis("Fe", "cr"), "species not available: Fe\\(cr\\)")
  ina <- nrow(thermo$obigt) + 1
  expect_error(basis(ina), "species not available")
  expect_error(preset.basis(c("CN")), "is not a keyword")
  # after all that, the basis should still be undefined
  expect_is(basis(), "NULL")
})

test_that("invalid basis modification requests cause an error", {
  basis(delete=TRUE)
  expect_error(mod.basis("CH4", "gas"), "basis is not defined")
  b <- basis("CHNOS+")
  expect_error(mod.basis("CH4", "gas"), "is not a formula of one of the basis species")
  iCH4 <- info("CH4")
  expect_error(mod.basis(iCH4, "gas"), "is not a species index of one of the basis species")
  expect_error(mod.basis("CO2", "PPM"), "the elements .* in buffer .* are not in the basis")
  expect_error(mod.basis("CO2", "liq"), "state .* not found")
  # after all that, the basis should be unchanged
  expect_equal(basis(), b)
})

test_that("changing state maintains species name", {
  b1 <- basis(c("Al2O3", "quartz", "oxygen"))
  b2 <- basis("SiO2", "cr2")
  # we went from quartz cr1 to cr2, which is the next row in the database
  expect_equal(sum(b2$ispecies - b1$ispecies), 1)
  expect_error(basis("SiO2", "cr3"), "state or buffer 'cr3' not found for quartz")
})
