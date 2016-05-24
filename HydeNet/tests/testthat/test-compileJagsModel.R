context("compileJagsModel")

data(PE, package="HydeNet")
Net <- HydeNetwork(~ wells +
                     pe | wells +
                     d.dimer | pregnant*pe +
                     angio | pe +
                     treat | d.dimer*angio +
                     death | pe*treat,
                   data = PE)

test_that("compileJagsModel returns an object of class 'compiledHydeNetwork'",
{
  compiledNet <- compileJagsModel(Net, n.chains=5)
  expect_equal(class(compiledNet),
               "compiledHydeNetwork")
})
