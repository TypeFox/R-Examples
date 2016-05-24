AICweights <- function(models){
  # models: vector of characters with the model names 
  nmodels <- length(models)
  AICc_values <- numeric(nmodels)
  for(i in 1:nmodels) AICc_values[i] <- AICc(eval(parse(text=models[i])))
  deltaAIC <- AICc_values - min(AICc_values)
  expdeltaAIC <- exp(-0.5*deltaAIC)
  AICweights <- expdeltaAIC/(sum(expdeltaAIC))
  names(AICweights) <- models
  return(AICweights)
}
