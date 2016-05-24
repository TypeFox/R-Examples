effectFactormod11 <-
function(object, level, nboot)
 {
  # Difference effect (needs bootstrap):
  auxd <- effectFactorBoot(object = object, level = level, nboot = nboot)
  effectDiff <- auxd$effect[, 1:3]
  if (!is.matrix(effectDiff)) effectDiff <- t(as.matrix(effectDiff))
  # Percent effect (doen't needs bootstrap):
  auxp <- effectFactorExactTrans(object = object, level = level)
  effectPerc <- 100 * (exp(auxp$effect) - 1)
  effect <- cbind(effectDiff, effectPerc)
  colnames(effect) <- rep(c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = "")), 2)
  colnames(effect)[c(1, 4)] <- paste(colnames(effect)[1], c("Diff", "Percent"), sep = "")  
  rownames(effect) <- auxd$rownameseffect
  info <- paste("Adjusted change in the mean of the response variable when the explanatory\nvariable changes from its reference level, '", auxd$Xbasal, "', to an alternative level\n(confidence interval for the difference change based on\n", nboot, " bootstrap samples)", sep = "")
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
