# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

rcr = function (formula, participants, dataframe, ...){
  if (length (participants) != nrow(dataframe)) return (cat("Error: Dataframe rows and participant vector length do not equal.\n\n"))
  parts = as.factor(levels (as.factor (participants)))
  nparts = length (parts)
  
  coefficients = matrix (0, nparts, ncol (model.matrix (formula, data = dataframe)))
  varExp = NULL  
  
  for (i in 1:nparts){
    temp = dataframe[participants == parts[i],]
    mod =  glm (formula, data = temp, ...)
    varExp = c(varExp, (mod$null.deviance - mod$deviance) / mod$null.deviance)
    coefficients[i,] = mod$coefficients
  }
  
  factor.names = attr(model.frame (formula, data = dataframe), "terms") 
  factor.names = c('(Intercept)', attr (factor.names, 'term.labels'))

  factors = attr(model.matrix (formula, data = dataframe), "assign")
  coefficients = data.frame (coefficients)
  coefficient.names = names(mod$coefficients)
  colnames(coefficients) = coefficient.names 
  coefficient.means = as.numeric (colMeans (coefficients))
  names (coefficient.means) = coefficient.names
  colnames(coefficients) = names(mod$coefficients)

  output =  (list (formula = formula, call = match.call(), participants = parts, 
  factors = as.factor(factors), factor.names = factor.names, coefficients = coefficients, 
  coefficient.means = coefficient.means, coefficient.names = coefficient.names, varExp = varExp))
  
  class (output) = 'rcr'
  output
}
