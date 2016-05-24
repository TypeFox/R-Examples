effectFactormod2 <-
function(object, level, nboot)
 {
  # Difference effect (doesn't needs bootstrap):
  auxd <- effectFactorExactTrans(object = object, level = level)
  effectDiff <- auxd$effect
  # Percent effect (needs bootstrap):
  auxp <- effectFactorBoot(object = object, level = level, nboot = nboot)
  effectPerc <- auxp$effect[, 4:6]
  if (!is.matrix(effectPerc)) effectPerc <- t(as.matrix(effectPerc))
  effect <- cbind(effectDiff, effectPerc)
  colnames(effect) <- rep(c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = "")), 2)
  colnames(effect)[c(1, 4)] <- paste(colnames(effect)[1], c("Diff", "Percent"), sep = "")  
  rownames(effect) <- auxd$rownameseffect
  info <- paste("Adjusted change in the mean of the response variable when the explanatory\nvariable changes from its reference level, '", auxd$Xbasal, "', to an alternative level\n(confidence interval for the percent change based on ", nboot, " bootstrap samples)", sep = "")
  if (any(MY(object)$M$Estimate < 0))
   info <- paste("WARNING: percent scale for effects not suitable because of negative values\nfor the adjusted mean.\n\n", info, sep = "")      
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
