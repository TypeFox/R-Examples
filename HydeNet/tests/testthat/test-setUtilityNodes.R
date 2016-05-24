context("setUtilityNodes")

test_that("setUtiltyNodes",
{
  expect_that(
    Net <- HydeNetwork(~ wells +
                       pe | wells +
                       d.dimer | pregnant*pe +
                       angio | pe +
                       treat | d.dimer*angio +
                       death | pe*treat,
                     data = PE) %>%
    setUtilityNodes(treat, angio),
    not(throws_error()))
})