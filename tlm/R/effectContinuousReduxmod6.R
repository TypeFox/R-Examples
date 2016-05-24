effectContinuousReduxmod6 <-
function(object, r, IQReffect, effectsOutRange, level)
 {
  # Percent effect (doesn't needs bootstrap):
  betaCI <- betaCIContinuousExactTrans(object = object, level = level)
  q <- 1 + r / 100
  effect <- 100 * (q^betaCI - 1)
  effect <- data.frame(r = r, effect)
  names(effect)[-1] <- c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = ""))
  if (IQReffect)
   {
   	info <- paste("Adjusted percent change in the geometric mean of the response variable\nfor an 'r'% change in the explanatory variable equivalent to the\ninterquartile ratio")
    } else {
    info <- paste("Adjusted percent change in the geometric mean of the response variable\nfor an 'r' = ", r, "% change in the explanatory variable", sep = "")
   }
  if (any(MY(object)$M$Estimate < 0))
   info <- paste("WARNING: percent scale for effects not suitable because of negative values\nfor the adjusted mean.\n\n", info, sep = "")
  if (!IQReffect & effectsOutRange)
   info <- paste("WARNING: computing effects out of the observed range of X.\n\n", info, sep = "")          
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
