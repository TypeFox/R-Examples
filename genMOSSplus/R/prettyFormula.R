prettyFormula <-
function (model, varNames) {

  model <- na.omit(model)
  formula <- paste("[", paste(varNames[model], collapse= ", "), "]", sep = "")

  return(formula)

}
