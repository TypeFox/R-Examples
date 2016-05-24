makeModelHierarchical <-
function (model, tools) {

  hierarchicalModel <- model

  for (i in 1:tools$nVarSets) {
    if (model[i] == 1)
      hierarchicalModel[getSubsets(i, tools)] <- 1
  }

  return(hierarchicalModel)

}
