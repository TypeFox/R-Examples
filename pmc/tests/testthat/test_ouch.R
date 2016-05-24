test_that("we can use ouch-based functions", {
  library("ouch")
  data(anoles)
  tree <- with(anoles, ouchtree(node, ancestor, time / max(time), species))

  ou3v4 <- pmc(tree, log(anoles["size"]), modelA = "hansen", modelB = "hansen", 
             optionsA = list(regimes = anoles["OU.LP"]), 
             optionsB = list(regimes = anoles["OU.4"]),
             nboot = 100, sqrt.alpha = 1, sigma = 1, mc.cores = 1)
  expect_is(ou3v4, "pmc")
})
