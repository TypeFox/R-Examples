getM <-
function(object, x, untransform, level)
 {  
  L <- getL(object = object, x = x)
  mod <- object$model 
  mf <- model.frame(mod)
  mt <- attr(mf, "terms")
  Xclass <- attr(mt, "dataClasses")[2]
  family <- family(mod)$family
  prob <- (1 + level) / 2
  if (any(family == c("binomial", "poisson")))
   z <- qnorm(prob) else  z <- qt(prob, df = mod$df.residual)
  L <- L["mean"] + c(0, -1, 1) * z * L["sd"]
  if (untransform)
   {
    M <- switch(family, gaussian = L, binomial = 1 / (1 + exp(-L)), poisson = exp(L))
    ypow <- object$ypow
    if (family == "gaussian" && ypow != 1)
     {
      Mbij <- bijectivityback(M, power = ypow)
      if (!Mbij)
         stop("non bijectivity of the computed expected values")
      M <- powerUntransform(M, ypow) 
     }
    } else {
    M <- L
   }
  names(M) <- c("Estimate", paste(c("lower", "upper"), 100 * level,  "%", sep = ""))
  if (Xclass != "factor")
   {
   	if (untransform)
     {
   	  M <- c(x, M)
   	  names(M)[1] <- "x"
   	  } else {
   	  xt <- powerTransform(x, object$xpow)
   	  M <- c(xt, M)
   	  names(M)[1] <- "xtrans"
   	  names(M)[2] <- paste(names(M)[2], "trans", sep = "")
     } } else {
   	Xlevels <- levels(mf[, 2])
   	xval <- which(Xlevels == x)
    M <- c(xval, M)
   	names(M)[1] <- "xlevel"
    if (!untransform)
     names(M)[2] <- paste(names(M)[2], "trans", sep = "")
   }	
  auxlo <- M[paste("lower", 100 * level,  "%", sep = "")]
  auxup <- M[paste("upper", 100 * level,  "%", sep = "")]
  if (any(auxlo > auxup))
   {
    M[paste("lower", 100 * level,  "%", sep = "")] <- auxup
    M[paste("upper", 100 * level,  "%", sep = "")] <- auxlo
   }
  M
 }
