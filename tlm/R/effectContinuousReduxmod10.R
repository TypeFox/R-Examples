effectContinuousReduxmod10 <-
function(object, c, IQReffect, effectsOutRange, level)
 {
  # Percent effect (doesn't needs bootstrap):
  betaCI <- betaCIContinuousExactTrans(object = object, level = level)
  effect <- 100 * (exp(betaCI * c) - 1)
  effect <- data.frame(c = c, effect)
  names(effect)[-1] <- c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = ""))
  if (IQReffect)
   {
   	info <- paste("Adjusted percent change in the mean of the response variable for a\n 'c' units additive change in the explanatory variable equivalent to\nthe interquartile range")
   	} else {
    info <- paste("Adjusted percent change in the mean of the response variable for a\n'c' =", c, "units additive change in the explanatory variable")
   }
  if (!IQReffect & effectsOutRange)
   info <- paste("WARNING: computing effects out of the observed range of X.\n\n", info, sep = "")          
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
