effectInfomod2 <-
function(object)
 {
  mf <- model.frame(object$model)
  aux <- summary(object$model)$coefficients
  Xlevels <- levels(mf[, 2])
  nlevels <- length(Xlevels)
  beta <- aux[2:nlevels, ]
  Xincrease <- paste("changing X from its reference, '", Xlevels[1], "', to the alternative level", sep = "")
  effecttype <- "additive change in the mean of Y"
  effectsize <- "beta units of Y"
  furtherinfo <- "\nFurther details can be obtained using effect() and providing the level for the\nconfidence interval, 'level'." 
  res <- list(beta = beta, Xincrease = Xincrease, effecttype = effecttype, effectsize = effectsize)
  res 	
 }
