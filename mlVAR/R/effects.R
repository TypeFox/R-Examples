

## fixedEffects:
fixedEffects <- function(object){
  # Nodes:
  Nodes <- rownames(object$fixedEffects)
  # Predictors:
  Predictor = colnames(object$fixedEffects)
  
  # Data frame of fixed effects:
  df <- data.frame(
    Response = Nodes[row(object$fixedEffects)],
    Predictor = Predictor[col(object$fixedEffects)],
    effect = c(object$fixedEffects),
    se = c(object$se.fixedEffects),
    p = c(object$pvals)
  )
  
  return(df)
}


# Random effects:
randomEffects<- function(object){
  # Nodes:
  Nodes <- rownames(object$randomEffectsVariance)
  # Predictors:
  Predictor = colnames(object$randomEffectsVariance)
  
  # Data frame of fixed effects:
  df <- data.frame(
    Response = Nodes[row(object$randomEffectsVariance)],
    Predictor = Predictor[col(object$randomEffectsVariance)],
    variance = c(object$randomEffectsVariance)
  )
  
  return(df)
}