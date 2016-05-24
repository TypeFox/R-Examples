effectFactormod8 <-
function(object, level, nboot)
 {
  # OR effect (doesn't needs bootstrap):
  aux <- effectFactorExactTrans(object = object, level = level)
  effect <- exp(aux$effect)
  info <- paste("Adjusted odds ratio of the response variable when the explanatory\nvariable changes from its reference level, '", aux$Xbasal, "', to an alternative level", sep = "")
  info <- paste(info, ":\n\n", sep = "")
  colnames(effect) <- c("OR", paste(c("lower", "upper"), 100 * level, "%", sep = ""))
  rownames(effect) <- aux$rownameseffect
  res <- list(effect = effect, info = info)
  res
 }
