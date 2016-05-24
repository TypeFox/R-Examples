findHierNeighbours <-
function (model, generators, dualGenerators, tools) {

  currentAddDel <- generators + dualGenerators
  nCurrentAddDel <- sum(currentAddDel == 1 & tools$lenVarSets >=2)
  neighbourList <- matrix(nrow = nCurrentAddDel, ncol = tools$nVarSets) 
  k <- 1

  for (i in 1:tools$nVarSets) {
    if (currentAddDel[i] == 1 && tools$lenVarSets[i] >= 2) {
      model[i] <- 1 - model[i]
      neighbourList[k,] <- model
      model[i] <- 1 - model[i] 
      k <- k + 1
    }
  } 

  return(neighbourList)

}
