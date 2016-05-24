effectContinuousmod9 <-
function(object, x1x2, level, nboot, IQReffect, effectsOutRange)
 {
  # OR effect (doesn't needs bootstrap):
  betaCI <- betaCIContinuousExactTrans(object = object, level = level)
  effect <- t(apply(x1x2, 1, FUN = function(x) (x[2] / x[1]) ^ betaCI))
  effect <- cbind(x1x2, effect)
  names(effect)[3:5] <- c("OR", paste(c("lower", "upper"), 100 * level, "%", sep = ""))
  if (IQReffect)
   {
   	info <- paste("Adjusted odds ratio of the response variable when the explanatory\nvariable changes from the 1st to the 3rd quartile", sep = "")
   	} else {
    info <- paste("Adjusted odds ratio of the response variable when the explanatory\nvariable changes from x1 to x2", sep = "")  		
   	}
  if (effectsOutRange)
   info <- paste("WARNING: computing effects out of the observed range of X.\n\n", info, sep = "")          
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
