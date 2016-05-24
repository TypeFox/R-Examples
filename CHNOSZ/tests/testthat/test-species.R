context("species")

# clear out any previous basis definition or database alterations
suppressMessages(data(thermo))

test_that("species not contained by basis cause errors", {
  expect_error(species("H2O"), "basis species are not defined")
  expect_error(species.basis("H2O"), "basis species are not defined")
  basis("CHNOS")
  expect_error(species.basis("U"), "element\\(s\\) not in the basis\\: U")
  expect_error(species("fayalite"), "element\\(s\\) not in the basis\\: Fe Si")
})

test_that("for one or more species, species.basis() keeps track of zeroes and puts elements in order of thermo$basis", {
  basis("CHNOS")
  test0 <- count.elements("OHN0")
  test1 <- count.elements("HN0O")
  expect_equal(species.basis(test0), species.basis(test1))
  # we can send multiple species to species.basis() but the argument has to be constructed correctly
  expect_equal(unique(as.numeric(species.basis(makeup(c("C", "CCN"))))), 0)
  expect_equal(species.basis(makeup(c("C", "CCN"), count.zero=TRUE))[2, , drop=FALSE], species.basis(makeup("CCN")))
})

test_that("deleting nonexistent species causes error or warning", {
  expect_error(species("CO2", delete=TRUE), "nonexistent species definition")
  species("H2O")
  expect_warning(species("CO2", delete=TRUE), "not present, so can not be deleted")
  expect_is(species("water", delete=TRUE), "NULL")
  # we should also get NULL if *all* species are deleted
  species("H2O")
  expect_is(species(delete=TRUE), "NULL")
})

test_that("non-available species cause error, and species can be added or modified", {
  basis("CHNOS")
  expect_error(species("wate"), "species not available")
  # add CO2, aq
  sdef <- species("CO2")
  # we can't add the same species twice
  expect_equal(nrow(species("CO2")), 1)
  # change it to gas
  expect_equal(species(1, "gas")$state, "gas")
  # change its log fugacity to -5
  expect_equal(species(1, -5)$logact, -5)
  # add CO2, aq
  expect_equal(nrow(species("CO2")), 2)
  # add alanine by index in thermo$obigt
  expect_equal(nrow(species(info("alanine"))), 3)
  # if we just use an index, get only that species
  expect_equal(species(3)$name, "alanine")
  # we can add a species with the same name but different state
  expect_equal(nrow(species("alanine", "cr")), 4)
  # we can modify the logact of a named species (only first match)
  expect_equal(species("alanine", -7)$logact[3], -7)
})

test_that("index_return provides indices for touched species", {
  basis("CHNOS")
  expect_equal(species("CO2", index.return=TRUE), 1)
  # here it's "touched" (but not added or modified)
  expect_equal(species("CO2", index.return=TRUE), 1)
  expect_equal(species(c("H2O", "NH3"), index.return=TRUE), c(2, 3))
  expect_equal(species(2, "gas", index.return=TRUE), 2)
})
