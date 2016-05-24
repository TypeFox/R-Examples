context("HydeNetwork")

data(PE, package='HydeNet')
Net <- HydeNetwork(~ wells
                       + pe | wells
                       + d.dimer | pregnant*pe
                       + angio | pe
                       + treat | d.dimer*angio
                       + death | pe*treat,
                       data = PE)

test_that("HydeNetwork.formula returns expected attributes",
{
  expect_equal(names(Net),
               c("nodes", "parents", "nodeType", "nodeFormula", "nodeFitter",
                 "nodeFitterArgs", "nodeParams", "fromData", "nodeData", 
                 "factorLevels",
                 "nodeModel", "nodeDecision", "nodePolicyValues", 
                 "nodeUtility", "dag", 
                 "data", "network_formula"))
})

test_that("HydeNetwork.formula assigns correct node types",
{
  expect_equal(Net$nodeType,
               list(wells = "dnorm", 
                    pe = "dbern",
                    d.dimer = "dnorm",
                    pregnant = "dcat",
                    angio = "dcat",
                    treat = "dbern",
                    death = "dcat"))
})

test_that("HydeNetwork.list returns expected attributes",
{
  g1 <- lm(wells ~ 1, data=PE)
  g2 <- glm(pe ~ wells, data=PE, family="binomial")
  g3 <- lm(d.dimer ~ pe + pregnant, data=PE)
  g4 <- xtabs(~ pregnant, data=PE)
  g5 <- cpt(angio ~ pe, data=PE)
  g6 <- glm(treat ~ d.dimer + angio, data=PE, family="binomial")
  g7 <- cpt(death ~ pe + treat, data=PE)
  
  bagOfModels <- list(g1,g2,g3,g4,g5,g6,g7)
  
  bagNet <- HydeNetwork(bagOfModels)
  
  expect_equal(names(bagNet),
               c("nodes", "parents", "nodeType", "nodeFormula", "nodeFitter",
                 "nodeFitterArgs", "nodeParams", "fromData", "nodeData",
                 "factorLevels", 
                 "nodeModel", "nodeDecision", "nodePolicyValues", "nodeUtility", "dag", 
                 "network_formula"))
})