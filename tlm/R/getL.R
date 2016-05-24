getL <-
function(object, x)
 {  
  mod <- object$model
  mf <- model.frame(mod)
  mt <- attr(mf, "terms")
  dessignMatrix <- model.matrix(mt, data = mf)
  X <- colMeans(dessignMatrix)
  Xclass <- attr(mt, "dataClasses")[2]
  if (Xclass != "factor")
   {
   	xt <- powerTransform(x, object$xpow)
    X[2] <- xt } else {
    xt <- x
    Xlevels <- levels(mf[, 2])
    whichlevel <- (levels(mf[, 2]) %in% xt)[-1]
    X[2:length(Xlevels)] <- whichlevel
   }
  aux <- summary(mod)$coefficients
  betas <- aux[, "Estimate"]
  Vbetas <- vcov(mod)	
  muL <- as.numeric(betas %*% X)
  varL <- as.numeric(t(X) %*% Vbetas %*% X)
  L <- c(muL, sqrt(varL))
  names(L) <- c("mean", "sd")
  L
 }
