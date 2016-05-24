context("compileDecisionModel")

Net <- HydeNetwork(~ wells +
                     pe | wells +
                     d.dimer | pregnant*pe +
                     angio | pe +
                     treat | d.dimer*angio +
                     death | pe*treat,
                   data = PE) %>%
  setDecisionNodes(angio, treat)

test_that("compileDecisionModel",
{
  expect_that(compileDecisionModel(Net),
              not(throws_error()))
})
