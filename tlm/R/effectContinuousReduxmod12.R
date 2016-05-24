effectContinuousReduxmod12 <-
function(object, r, IQReffect, effectsOutRange, level)
 {
  # Percent effect (doesn't needs bootstrap):
  betaCI <- betaCIContinuousExactTrans(object = object, level = level)
  q <- 1 + r / 100
  effect <- 100 * (q ^ betaCI - 1)
  effect <- data.frame(r = r, effect)
  names(effect)[-1] <- c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = ""))
  if (IQReffect)
   {
   	info <- paste("Adjusted percent change in the mean of the response variable for an\n'r'% change in the explanatory variable equivalent to the interquartile\nratio", sep = "")
   	} else {
   	info <- paste("Adjusted percent change in the mean of the response variable for an 'r' = ", r, "%\nchange in the explanatory variable", sep = "")
   }
  if (!IQReffect & effectsOutRange)
   info <- paste("WARNING: computing effects out of the observed range of X.\n\n", info, sep = "")          
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
