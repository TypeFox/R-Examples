context("HydePosterior")

data(PE, package="HydeNet")
Net <- HydeNetwork(~ wells +
                     pe | wells +
                     d.dimer | pregnant*pe +
                     angio | pe +
                     treat | d.dimer*angio +
                     death | pe*treat,
                   data = PE)
compiledNet <- compileJagsModel(Net, n.chains=5)

test_that("Unbound HydePosterior returns object of class HydePosterior",
{
  Posterior <- HydePosterior(compiledNet,
                             variable.names = c("d.dimer", "death"),
                             n.iter = 1000, 
                             bind = FALSE)  
  expect_equal(class(Posterior),
               "HydePosterior")
})

test_that("Bound HydePosterior returns object of class data.frame",
{
  Posterior <- HydePosterior(compiledNet,
                             variable.names = c("d.dimer", "death"),
                             n.iter = 1000)  
  expect_equal(class(Posterior),
               "data.frame")
})

test_that("Unbound HydePosterior print method",
{
  Posterior <- HydePosterior(compiledNet,
                             variable.names = c("d.dimer", "death"),
                             n.iter = 1000, 
                             bind = FALSE)  
  expect_that(print(Posterior), 
              not(throws_error()))
})

test_that("bindPosterior returns relabeled data",
{
  Posterior <- HydePosterior(compiledNet,
                             variable.names = c("d.dimer", "death"),
                             n.iter = 1000, 
                             bind = FALSE) 
  Bound <- bindPosterior(Posterior)
  expect_equal(class(Bound$death),
               "factor")
})

test_that("bindPosterior returns relabeled data",
{
  Posterior <- HydePosterior(compiledNet,
                             variable.names = c("d.dimer", "death"),
                             n.iter = 1000, 
                             bind = FALSE) 
  Bound <- bindPosterior(Posterior, relabel_factor = FALSE)
  expect_equal(class(Bound$death),
               "numeric")
})