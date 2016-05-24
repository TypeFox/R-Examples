context("factorFormula")

test_that("factorFormula",
{
  Net <- HydeNetwork(~ wells +
                       pe | wells +
                       d.dimer | pregnant*pe +
                       angio | pe +
                       treat | d.dimer*angio +
                       death | pe*treat,
                     data = PE)
  expect_equal(factorFormula(death ~ ilogit((treat == "No") + (angio == "Positive")),
                             Net),
               death ~ ilogit((treat == 0) + (angio == 2)))
})
  