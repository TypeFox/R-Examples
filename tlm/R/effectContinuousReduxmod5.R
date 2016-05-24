effectContinuousReduxmod5 <-
function(object, r, IQReffect, effectsOutRange, level)
 {
  # Difference effect (doesn't needs bootstrap):
  betaCI <- betaCIContinuousExactTrans(object = object, level = level)
  q <- 1 + r / 100
  effect <- betaCI * log(q)
  effect <- data.frame(r = r, effect)
  names(effect)[-1] <- c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = ""))
  if (IQReffect)
   {
   	info <- paste("Adjusted additive change in the mean of the response variable for an\n'r'% change in the explanatory variable equivalent to the interquartile\nratio")
    } else {
    info <- paste("Adjusted additive change in the mean of the response variable for an\n'r' = ", r, "% change in the explanatory variable", sep = "")
   }
  if (!IQReffect & effectsOutRange)
   info <- paste("WARNING: computing effects out of the observed range of X.\n\n", info, sep = "")          
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
