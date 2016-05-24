effectFactor <-
function(object, level, nboot, modeltype)
 { 
  # No/Yes bootstrap cases (depending on the effect scale):
  if(any(modeltype %in% c(2, 4, 8, 11)))
   {
   	eval(parse(text = paste("res <- effectFactormod", modeltype, "(object = object, level = level, nboot = nboot)", sep = "")))
   	} else {
    # Bootstrap cases:
    aux <- effectFactorBoot(object = object, level = level, nboot = nboot)
    effect <- aux$effect
    info <- paste("Adjusted change in the median of the response variable when the explanatory\nvariable changes from its reference level, '", aux$Xbasal, "', to an alternative level\n(based on ", nboot, " bootstrap samples)", sep = "")    
    if (any(MY(object)$M$Estimate < 0))
     info <- paste("WARNING: percent scale for effects not suitable because of negative values\nfor the adjusted median.\n\n", info, sep = "")    
    info <- paste(info, ":\n\n", sep = "")
    colnames(effect) <- rep(c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = "")), 2)
    colnames(effect)[c(1, 4)] <- paste(colnames(effect)[1], c("Diff", "Percent"), sep = "")  
    rownames(effect) <- aux$rownameseffect
    res <- list(effect = effect, info = info)
   }
  res
 }
