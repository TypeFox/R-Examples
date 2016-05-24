effectFactorExactTrans <-
function(object, level)
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
  effect <- coefs[2:nlevels, 1]
  CI <- confint(mod, level = level)[2:nlevels, ]
  if (nlevels == 2)
   {
    effect <- matrix(c(effect, CI), nrow = 1, ncol = 3, byrow = TRUE)
    names(effect)[1] <- "Estimate"
    } else {
    effect <- cbind(effect, CI)
   }
  res <- list(effect = effect, Xbasal = Xlevels[1], rownameseffect = rownameseffect)
  res
 }
