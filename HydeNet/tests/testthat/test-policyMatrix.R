context("policyMatrix")

Net <- HydeNetwork(~ wells +
                     pe | wells +
                     d.dimer | pregnant*pe +
                     angio | pe +
                     treat | d.dimer*angio +
                     death | pe*treat,
                   data = PE) %>%
  setDecisionNodes(angio, treat)

test_that("Get the default policy matrix",
{
  expect_equal(policyMatrix(Net),
               expand.grid(angio = 1:2, treat = 0:1))
})

test_that("Get a customized policy matrix",
{
  expect_equal(policyMatrix(Net, angio = 1, treat = 0:1),
               expand.grid(angio = 1, treat = 0:1))
})

test_that("Get a customized policy matrix with a continuous variable",
{
  expect_equal(policyMatrix(Net, angio = 1, treat = 0:1, d.dimer = 3),
               expand.grid(angio = 1, treat = 0:1, d.dimer = 3))
})

test_that("Cast error when using a variable not in the network",
{
  expect_error(policyMatrix(Net, angio = 1, treat = 0:1, x = 2))
})