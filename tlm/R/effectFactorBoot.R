effectFactorBoot <-
function(object, level, nboot)
 {
  mod <- object$model
  mf <- model.frame(mod)
  coefs <- summary(mod)$coefficients
  Xlevels <- levels(mf[, 2])
  nlevels <- length(Xlevels)
  x1 <- rep(Xlevels[1], nlevels - 1)
  x2 <- Xlevels[-1]
  x2lab <- paste(x2, ":", sep = "")
  rownameseffect <- paste(x1, x2, sep = " -> ")
  aux <- data.frame(x1 = x1, x2 = x2)
  effect <- t(apply(aux, 1, FUN = function(x) geteffectx1x2(object = object, x1 = x[1], x2 = x[2], level = level, nboot = nboot)))
   if (!is.matrix(effect))
  # effect <- matrix(effect, nrow = 1, ncol = 3, byrow = TRUE)
    effect <- matrix(effect, nrow = 1, ncol = 6, byrow = TRUE)
  res <- list(effect = effect, Xbasal = Xlevels[1], rownameseffect = rownameseffect)
  res
 }
