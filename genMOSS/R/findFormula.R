findFormula <-
function (generators, tools) {

  if (generators[1] == 1) {
    formula <- paste ("freq", " ~ ", 1, sep = "")
    return(formula)
  }

  formula <- paste ("freq", " ~ ", sep = "")

  for(i in tools$nVarSets:1) {
    if(generators[i] == 1) {
      first <- 1
      for(j in 1:tools$n) {
        if(tools$varSets[i,j] == 1) { 
          if (first == 1) {
            first <- 0  
            formula <- paste (formula,  tools$varNames[j], sep = "")
          }
          else 
            formula <- paste (formula, "*", tools$varNames[j], sep = "")
        }
      }

      formula <- paste (formula, " + ", sep = "")    
    }
  } 
  
  formula <- substr (formula, 1,nchar(formula)-3)
  return(formula)    

}
