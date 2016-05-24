findDualGenerators <-
function (model, tools) {

  dualGenerators <- 1 - model
  
  for (i in 1:tools$nVarSets) {
    if (model[i] == 0) {
      index <- na.omit(tools$upLinks[i,])      
      dualGenerators[index] <- 0      
    }
  }

  return(dualGenerators)

}
