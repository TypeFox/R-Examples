modelType <-
function(object)
 {
  mod <- object$model
  family <- family(mod)$family
  xpow <- object$xpow
  ypow <- object$ypow
  mf <- model.frame(mod)
  mt <- attr(mf, "terms")
  Xclass <- attr(mt, "dataClasses")[2]
  nmodels <- 12
  cond <- rep(NA, nmodels)
  cond[1] <- family == "gaussian" && Xclass != "factor" && ypow == 1 && xpow == 1
  cond[2] <- family == "gaussian" && Xclass == "factor" && ypow == 1
  cond[3] <- family == "gaussian" && Xclass != "factor" && ypow == 0 && xpow == 1
  cond[4] <- family == "gaussian" && Xclass == "factor" && ypow == 0
  cond[5] <- family == "gaussian" && Xclass != "factor" && ypow == 1 && xpow == 0
  cond[6] <- family == "gaussian" && Xclass != "factor" && ypow == 0 && xpow == 0
  cond[7] <- family == "binomial" && Xclass != "factor" && xpow == 1
  cond[8] <- family == "binomial" && Xclass == "factor" && xpow == 1
  cond[9] <- family == "binomial" && Xclass != "factor" && xpow == 0
  cond[10] <- family == "poisson" && Xclass != "factor" && xpow == 1
  cond[11] <- family == "poisson" && Xclass == "factor" && xpow == 1
  cond[12] <- family == "poisson" && Xclass != "factor" && xpow == 0
  modeltype <- as.numeric(1:nmodels %*% cond)
  modeltype
 }
