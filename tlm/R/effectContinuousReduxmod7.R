effectContinuousReduxmod7 <-
function(object, c, IQReffect, effectsOutRange, level)
 {
  # OR effect (doesn't needs bootstrap):
  betaCI <- betaCIContinuousExactTrans(object = object, level = level)
  effect <- exp(c * betaCI)
  effect <- data.frame(c = c, effect)
  names(effect)[-1] <- c("OR", paste(c("lower", "upper"), 100 * level, "%", sep = ""))
  if (IQReffect)
   {
    info <- paste("Adjusted odds ratio of the response variable for a 'c' units additive\nchange in the explanatory variable equivalent to the interquartile range")
    } else {
    info <- paste("Adjusted odds ratio of the response variable for a 'c' =", c, "\nunits additive change in the explanatory variable")
   }
  if (!IQReffect & effectsOutRange)
   info <- paste("WARNING: computing effects out of the observed range of X.\n\n", info, sep = "")          
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
