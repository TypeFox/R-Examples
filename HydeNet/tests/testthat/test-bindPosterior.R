context("bindPosterior")

data(PE, package="HydeNet")
Net <- HydeNetwork(~ wells +
                     pe | wells +
                     d.dimer | pregnant*pe +
                     angio | pe +
                     treat | d.dimer*angio +
                     death | pe*treat,
                   data = PE) %>%
  setDecisionNodes(treat, angio)
Post <- compileDecisionModel(Net) %>%
  HydePosterior(variable.names = c("wells", "treat", "death"),
                n.iter = 100,
                bind = FALSE) 

test_that("bindPosterior from Decision Model",
{
  expect_that(bindPosterior(Post), not(throws_error()))
})
