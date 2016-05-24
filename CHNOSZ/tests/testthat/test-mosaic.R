context("mosaic")

test_that("results are consistent with affinity()", {
  basis(c("CO2", "H2O", "NH3", "O2"), c(0, 0, 0, 0))
  species(c("alanine", "glycine"))
  a <- affinity()
  # this is a degenerate case because we only allow NH3 to swap for NH3, and CO2 for CO2;
  # however it still exercises the affinity scaling and summing code
  m1 <- mosaic("NH3", "CO2", blend=TRUE)
  # this failed before we divided by loga.tot to get _relative_ abundances of basis species in mosaic.R
  expect_equal(a$values, m1$A.species$values)
  # the next call failed when which.pmax(), called by diagram(), choked on a list of length one
  m2 <- mosaic("NH3", "CO2")
  expect_equal(a$values, m2$A.species$values)
})

test_that("blend=TRUE produces reasonable values", {
  # a more rigorous test than above. this was failing because loga.tot (actually, a.tot)
  # was computed incorrectly, by sum'ing an unlist'ed list (the affinities of basis species)
  # to produce a single value; corrected by using Reduce for addition of vectors/arrays in the list.
  # example adapted from ?mosaic
  basis(c("FeO", "SO4-2", "H2O", "H+", "e-"))
  basis("SO4-2", -6)
  basis("Eh", -0.15)
  species(c("hematite", "magnetite"))
  # the basis species we'll swap through
  bases <- c("SO4-2", "HSO4-", "HS-", "H2S")         
  # calculate affinities using the predominant basis species
  pH <- c(0, 14, 29)
  m1 <- mosaic(bases, pH=pH)
  m2 <- mosaic(bases, pH=pH, blend=TRUE)
  # these species have no S so the results should be the same
  expect_equal(m1$A.species$values, m2$A.species$values)
  species(c("pyrrhotite", "pyrite"))
  # now with S-bearing species ...
  m3 <- mosaic(bases, pH=pH)
  m4 <- mosaic(bases, pH=pH, blend=TRUE)
  # the results are different ...
  expect_equal(sapply(m3$A.species$values, "[", 13), sapply(m4$A.species$values, "[", 13), tol=1e-1)
  # but more similar at extreme pH values
  expect_equal(sapply(m3$A.species$values, "[", 1), sapply(m4$A.species$values, "[", 1), tol=1e-6)
  expect_equal(sapply(m3$A.species$values, "[", 29), sapply(m4$A.species$values, "[", 29), tol=1e-10)
})

# TODO: test that basis specifications can be exchanged between bases and bases2 without altering output
