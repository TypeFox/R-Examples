effectFactormod4 <-
function(object, level, nboot)
 {
  # Difference effect (needs bootstrap):
  auxd <- effectFactorBoot(object = object, level = level, nboot = nboot)
  effectDiff <- auxd$effect[, 1:3]
  if (!is.matrix(effectDiff)) effectDiff <- t(as.matrix(effectDiff))
  # Percent effect (doesn't needs bootstrap):
  auxp <- effectFactorExactTrans(object = object, level = level)
  effectPerc <- 100 * (exp(auxp$effect) - 1)
  effect <- cbind(effectDiff, effectPerc)
  colnames(effect) <- rep(c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = "")), 2)
  colnames(effect)[c(1, 4)] <- paste(colnames(effect)[1], c("Diff", "Percent"), sep = "")  
  rownames(effect) <- auxd$rownameseffect
  info <- paste("Adjusted change in the geometric mean of the response variable when\nthe explanatory variable changes from its reference level, '",auxd$Xbasal, "', to\nan alternative level (confidence interval for the difference based\non ", nboot, " bootstrap samples)", sep = "")
  if (any(MY(object)$M$Estimate < 0))
   info <- paste("WARNING: percent scale for effects not suitable because of negative values\nfor the adjusted mean.\n\n", info, sep = "")      
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
