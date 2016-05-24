prettyFormula <-
function (model, varNames) {

  model <- na.omit(model)
  
  formula <- paste("[", varNames[model[length(model)]], " | ", sep = "")
  formula <- paste(formula, paste(varNames[model[-length(model)]], collapse= ", "), "]", sep = "")
  #formula <- paste("[", paste(varNames[model], collapse= ", "), "]", sep = "")

  return(formula)

}
