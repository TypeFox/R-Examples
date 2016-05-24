effectContinuousReduxmod9 <-
function(object, r, IQReffect, effectsOutRange, level)
 {
  # OR effect (doesn't needs bootstrap):
  betaCI <- betaCIContinuousExactTrans(object = object, level = level)
  q <- 1 + r / 100
  effect <- q ^ betaCI
  effect <- data.frame(r = r, effect)
  names(effect)[-1] <- c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = ""))
  if (IQReffect)
   {
   	info <- paste("Adjusted odds ratio of the response variable for an 'r'% change in\nthe explanatory variable equivalent to the interquartile ratio", sep = "")
   	} else {
    info <- paste("Adjusted odds ratio of the response variable for an 'r' = ", r, "% change\nin the explanatory variable", sep = "")
   }
  if (!IQReffect & effectsOutRange)
   info <- paste("WARNING: computing effects out of the observed range of X.\n\n", info, sep = "")          
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
