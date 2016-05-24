effectContinuousReduxmod3 <-
function(object, c, IQReffect, effectsOutRange, level)
 {
   # Percent effect (doesn't needs bootstrap):
  betaCI <- betaCIContinuousExactTrans(object = object, level = level)
  effect <- 100 * (exp(betaCI * c) - 1)
  effect <- data.frame(c = c, effect)
  names(effect)[-1] <- c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = ""))
  if (IQReffect)
   {
  	info <- paste("Adjusted percent change in the geometric mean of the response variable\nfor a 'c' units additive change in the explanatory variable equivalent\nto the interquartile range")
    } else {
    info <- paste("Adjusted percent change in the geometric mean of the response variable\nfor a 'c' =", c, "units additive change in the explanatory variable")
   }
  if (any(MY(object)$M$Estimate < 0))
   info <- paste("WARNING: percent scale for effects not suitable because of negative values\nfor the adjusted mean.\n\n", info, sep = "")
  if (!IQReffect & effectsOutRange)
   info <- paste("WARNING: computing effects out of the observed range of X.\n\n", info, sep = "")          
  info <- paste(info, ":\n\n", sep = "")
  res <- list(effect = effect, info = info)
  res
 }
