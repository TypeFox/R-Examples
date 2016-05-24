findModel <-
function (formula, tools) {

  formula <- as.formula(formula)
  terms <- terms(formula)
  terms <- attributes(terms)$factors
  terms <- terms[-1,]
  terms <- t(terms)
  terms <- terms[,tools$varNames]  
  model <- rep(0,tools$nVarSets)  

  for (i in 1:dim(terms)[1])
    model[binToDec(terms[i,]) + 1] <- 1

  model[1] <- 1  

  return(model)  

}
