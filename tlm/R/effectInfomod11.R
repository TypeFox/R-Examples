effectInfomod11 <-
function(object)
 {
  mf <- model.frame(object$model)
  aux <- summary(object$model)$coefficients
  Xlevels <- levels(mf[, 2])
  nlevels <- length(Xlevels)
  beta <- aux[2:nlevels, ]
  Xincrease <- paste("changing X from its reference, '", Xlevels[1], "', to the alternative level", sep = "")
  effecttype <- "percent change in the mean of Y"
  effectsize <- "100 * [exp(beta) - 1]%"
  res <- list(beta = beta, Xincrease = Xincrease, effecttype = effecttype, effectsize = effectsize)
  res
 }
