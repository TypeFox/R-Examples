context("setNodeModels")

Net <- HydeNetwork(~ wells +
                     pe | wells +
                     d.dimer | pregnant*pe +
                     angio | pe +
                     treat | d.dimer*angio +
                     death | pe*treat)

g1 <- lm(wells ~ 1, data=PE)
g2 <- glm(pe ~ wells, data=PE, family="binomial")
g3 <- lm(d.dimer ~ pe + pregnant, data=PE)
g4 <- xtabs(~ pregnant, data=PE)
g5 <- cpt(angio ~ pe, data=PE)
g6 <- glm(treat ~ d.dimer + angio, data=PE, family="binomial")
g7 <- cpt(death ~ pe + treat, data=PE)

bagOfModels <- list(g1,g2,g3,g4,g5,g6,g7)

test_that("Returns a network",
{
  expect_that(setNodeModels(Net, g1,g2,g3,g4,g5,g6,g7),
              not(throws_error()))
})

test_that("Cast error when no models given",
{
  expect_error(setNodeModels(Net))
})

test_that("Cast error when not applying to a HydeNetwork",
{
  expect_error(setNodeModels(bagOfModels, g1))
})

test_that("Response variable is a node in the HydeNetwork",
{
  expect_error(setNodeModels(Net, lm(mpg ~ am, data = mtcars)))
})

test_that("Check that all regression variables are parents of the response",
{
  expect_error(setNodeModels(Net, lm(d.dimer ~ pe + pregnant + treat, data = PE)))
})
  