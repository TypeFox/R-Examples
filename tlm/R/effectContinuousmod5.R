effectContinuousmod5 <-
function(object, x1x2, level, nboot, IQReffect, effectsOutRange)
 {
  # Difference effect (doesn't needs bootstrap):
  betaCI <- betaCIContinuousExactTrans(object = object, level = level)
  effectDiff <- t(apply(x1x2, 1, FUN = function(x) betaCI * log(x[2] / x[1])))
  # Percent effect (needs bootstrap):
  auxp <- t(apply(x1x2, 1, FUN = function (x) geteffectx1x2(object = object, x1 = x[1], x2 = x[2], level = level, nboot = nboot)))
  effectPerc <- auxp[, 4:6]
  if (!is.matrix(effectPerc)) effectPerc <- t(as.matrix(effectPerc))
  effect <- cbind(x1x2, effectDiff, effectPerc)
  colnames(effect) <- c(paste("x", 1:2, sep = ""), rep(c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = "")), 2))
  colnames(effect)[c(3, 6)] <- paste(names(effect)[3], c("Diff", "Percent"), sep = "")
  if (IQReffect)
   {
   	info <- paste("Adjusted change in the mean of the response variable when the explanatory\nvariable changes from the 1st to the 3rd quartile (confidence interval\nfor the percent change based on ", nboot, " bootstrap samples)", sep = "")
   	} else {
    info <- paste("Adjusted change in the mean of the response variable when the explanatory\nvariable changes from x1 to x2 (confidence interval for the percent change\nbased on ", nboot, " bootstrap samples)", sep = "")  		
   	}
  if (any(MY(object)$M$Estimate < 0))
   info <- paste("WARNING: percent scale for effects not suitable because of negative values\nfor the adjusted mean.\n\n", info, sep = "")
  if (effectsOutRange)
   info <- paste("WARNING: computing effects out of the observed range of X.\n\n", info, sep = "")          
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
