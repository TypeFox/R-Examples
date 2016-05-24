context("expectedVariables")

test_that("expectedVariables",
{
  Net <- HydeNetwork(~ wells +
                       pe | wells +
                       d.dimer | pregnant*pe +
                       angio | pe +
                       treat | d.dimer*angio +
                       death | pe*treat,
                     data = PE) %>%
    setDecisionNodes(treat, angio)
  expect_equal(expectedVariables(Net, treat, TRUE),
               c("d.dimer", "angio"))
})