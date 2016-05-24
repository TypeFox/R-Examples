`fit.Models.glmFormula` <-
function(fmla, data.Vector, train.Frame, family = binomial, ...){
  glm(formula = fmla, data=train.Frame, family = family)
}

