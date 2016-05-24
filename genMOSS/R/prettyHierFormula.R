prettyHierFormula <-
function (generators, tools) {
 
  formula <- c()

  for(i in tools$nVarSets:1) {
    if(generators[i] == 1) {
      formula <- paste (formula, "[", sep = "")
      first <- 1
      for(j in 1:tools$n) {
        if(tools$varSets[i,j] == 1) { 
          if (first == 1) {
            first <- 0  
            formula <- paste (formula,  tools$varNames[j], sep = "")
          }
          else 
            formula <- paste (formula, ",", tools$varNames[j], sep = "")
        }
      }

      formula <- paste (formula, "]", sep = "")    
    }
  } 
  
  return(formula)

}
