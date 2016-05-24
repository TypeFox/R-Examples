findGenerators <-
function (model, tools) {

  generators <- model  

  for (i in 1:tools$nVarSets) {
    if (model[i] == 1) {
      index <- na.omit(tools$downLinks[i,])      
      generators[index] <- 0      
    }
  }
 
  return(generators)

}
