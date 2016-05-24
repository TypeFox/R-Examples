randomHierModel <-
function (p, tools) {

  model <- rbinom(tools$nVarSets, 1, p)
  model[tools$lenVarSets == 1] <- 1
  model <- makeModelHierarchical(model, tools)
  return(model)

}
