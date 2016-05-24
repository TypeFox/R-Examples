#I'm not thoroughly testing the output of the quantile, mantel, etc. methods, because those are maintained by other people
# - though perhaps a few worked examples? Might help to check my own code...
#Setup
require(pez)
require(testthat)
require(picante)
data(phylocom)
data <- comparative.comm(phylocom$phylo, phylocom$sample, traits=phylocom$traits, warn=FALSE)

context("fingerprint.regression")

test_that("quantile", {
  set.seed(123)
  basic.quantile <<- fingerprint.regression(data=data, eco.permute=100)
  expect_that(basic.quantile, is_a("fingerprint.regression"))
  expect_that(basic.quantile$eco.method, equals("quantile"))
  expect_that(basic.quantile$evo.method, equals("lambda"))
  expect_that(names(basic.quantile), equals(c("evo", "eco", "evo.method", "eco.method")))
  expect_that(basic.quantile$evo, is_equivalent_to(c(1, 1, 1e-06, 1)))
  set.seed(123)
  expect_that(basic.quantile$eco$raw, equals(eco.trait.regression(data, "taxa.labels", 100, "quantile", altogether=FALSE)))
  set.seed(123)
  expect_that(basic.quantile$evo, equals(phy.signal(data, method="lambda")))
})

test_that("lm", {
  set.seed(123)
  basic.lm <<- fingerprint.regression(data, eco.permute=100, eco.method="lm", eco.rnd="richness", evo.method="delta")
  expect_that(basic.lm, is_a("fingerprint.regression"))
  expect_that(basic.lm$eco.method, equals("lm"))
  expect_that(basic.lm$evo.method, equals("delta"))
  expect_that(names(basic.lm), equals(c("evo", "eco", "evo.method", "eco.method")))
  set.seed(123)
  expect_that(basic.lm$eco$raw, equals(eco.trait.regression(data, "richness", 100, "lm", altogether=FALSE)))
  set.seed(123)
  expect_that(basic.lm$evo, equals(phy.signal(data, method="delta")))
})

test_that("mantel", {
  set.seed(123)
  basic.mantel <<- fingerprint.regression(data, eco.permute=10, eco.method="mantel", eco.rnd="trialswap", evo.method="kappa")
  expect_that(basic.mantel, is_a("fingerprint.regression"))
  expect_that(basic.mantel$eco.method, equals("mantel"))
  expect_that(basic.mantel$evo.method, equals("kappa"))
  expect_that(names(basic.mantel), equals(c("evo", "eco", "evo.method", "eco.method")))
  set.seed(123)
  expect_that(basic.mantel$eco$raw, equals(eco.trait.regression(data, "trialswap", 10, "mantel", altogether=FALSE)))
  set.seed(123)
  expect_that(basic.mantel$evo, equals(phy.signal(data, method="kappa")))
})
