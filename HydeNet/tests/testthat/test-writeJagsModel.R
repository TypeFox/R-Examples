context("writeJagsModel")

craps <- HydeNetwork(~ d1 + d2 + diceSum | d1*d2
                     + firstRollOutcome | diceSum)  %>%
  setNode(d1, nodeType="dcat",
          pi = vectorProbs(p = rep(1/6,6), d1),
          validate = FALSE) %>%
  setNode(d2, nodeType="dcat",
          pi = vectorProbs(p = rep(1/6,6), d2),
          validate = FALSE) %>%
  setNode(diceSum, nodeType = "determ",
          define = fromFormula(),
          nodeFormula = diceSum ~ di1 + di2) %>%
  setNode(firstRollOutcome, nodeType = "determ",
          define = fromFormula(),
          nodeFormula = firstRollOutcome ~ 
            ifelse(diceSum < 4 | diceSum > 11, -1,
                   ifelse(diceSum == 7 | diceSum == 11, 1,0)))

test_that("writeJagsModel - determ",
{
  #* Because of how `writeJagsModel` processes the network object 
  #* it needs to be used within `writeNetworkModel`, which is 
  #* intended anyway since it is an unexported function
  expect_that(writeNetworkModel(craps, TRUE),
              not(throws_error()))
})

test_that("writeJagsModel - dcat with pi defined by user",
{
  expect_that(writeNetworkModel(craps, TRUE),
              not(throws_error()))
})

test_that("writeJagsModel - dpois",
{
  carNet <- HydeNetwork(~gear | mpg + am,
                        data = mtcars) %>%
  setNode(gear, nodeType = "dpois", nodeFitter = "glm",
            fitterArgs = list(family = poisson),
            lambda = fromData())
  expect_that(writeNetworkModel(carNet, TRUE),
              not(throws_error()))
})

