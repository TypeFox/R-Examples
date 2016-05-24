effectContinuousmod10 <-
function(object, x1x2, level, nboot, IQReffect, effectsOutRange)
 {
  # Difference effect (needs bootstrap):
  auxd <- t(apply(x1x2, 1, FUN = function (x) geteffectx1x2(object = object, x1 = x[1], x2 = x[2], level = level, nboot = nboot)))
  effectDiff <- auxd[, 1:3]
  if (!is.matrix(effectDiff)) effectDiff <- t(as.matrix(effectDiff))
  # Percent effect (doesn't needs bootstrap):
  betaCI <- betaCIContinuousExactTrans(object = object, level = level)
  effectPerc <- t(apply(x1x2, 1, FUN = function(x) 100 * (exp(betaCI * (x[2] - x[1])) - 1)))
  effect <- cbind(x1x2, effectDiff, effectPerc)
  colnames(effect) <- c(paste("x", 1:2, sep = ""), rep(c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = "")), 2))
  colnames(effect)[c(3, 6)] <- paste(names(effect)[3], c("Diff", "Percent"), sep = "")
  if (IQReffect)
   {
   	info <- paste("Adjusted change in the mean of the response variable when the explanatory\nvariable changes from the 1st to the 3rd quartile\n(confidence interval for the difference change based on\n", nboot, " bootstrap samples)", sep = "")
   	} else {
    info <- paste("Adjusted change in the mean of the response variable when the explanatory\nvariable changes from x1 to x2\n(confidence interval for the difference change based on\n", nboot, " bootstrap samples)", sep = "")  		
   	}
  if (effectsOutRange)
   info <- paste("WARNING: computing effects out of the observed range of X.\n\n", info, sep = "")          
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
