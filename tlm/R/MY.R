MY <-
function(object, x = NULL, npoints = 10, space = c("original", "transformed"), level = 0.95)
 { 
  if (!inherits(object, "tlm"))
     stop("argument 'object' must be of class 'tlm'")
    
  if(any(is.na(coef(object$model))))
     stop("MY is not available for models with any missing estimated coefficient")
  
  if (!is.null(x) && !inherits(x, "numeric") && !inherits(x, "integer"))
    stop("'x' must be a number or a numeric vector")
  
  ypow <- object$ypow
  xpow <- object$xpow
  
  if (!is.null(x) && xpow != 1)
    {
     xbij <- bijectivityforward(x, power = xpow)
      if (!xbij)
         stop("non bijectivity of the provided 'x'")
    }
  
  if ((!inherits(npoints, "numeric") && !inherits(npoints, "integer")) || npoints <= 0 || length(npoints) != 1)
    stop("'npoints' must be a positive integer")
     
  space <- match.arg(space) 
  
  if (!inherits(level, "numeric")  || level <= 0 || level >= 1 || length(level) != 1)
    stop("'level' must be a number in (0, 1)")
    
  if (space == "transformed" && family(object$model)$family == "gaussian" && xpow == 1 && object$ypow == 1)
    stop("there is no transformations in this model")
   
  mod <- object$model
  mf <- model.frame(mod)
  mt <- attr(mf, "terms")
  dessignMatrix <- model.matrix(mt, data = mf)
  Xclass <- attr(mt, "dataClasses")[2]
  cond1 <- Xclass != "factor"
  cond2 <- is.null(x) || !inherits(x, "numeric") && !inherits(x, "integer")
  cond3 <- is.null(npoints) || !inherits(npoints, "numeric") && !inherits(npoints, "integer") || length(npoints) != 1L
  if (cond1 && cond1 && cond3)
   stop("either a numeric value (or vector) 'x' or the number of points 'npoints' must be provided")
   
  if (Xclass == "factor")
   {
   	xvector <- levels(mf[, 2])
    } else {
    if (!is.null(x))
     xvector <- as.list(x) else
     {
      # Puntos equiespaciados en el espacio correspondiente:
      xtransmin <- min(dessignMatrix[, 2])	
      xtransmax <- max(dessignMatrix[, 2])
      if (space == "original")
       # secuencia de npoints puntos equidistantes en original
       xvector <- seq(powerUntransform(xtransmin, xpow), powerUntransform(xtransmax, xpow), length.out = npoints) else
       # secuencia de npoints puntos equidistantes en transformed y destransformar
       xvector <- powerUntransform(seq(xtransmin, xtransmax, length.out = npoints), xpow)
      }	
   }
  res <- as.data.frame(t(sapply(xvector, FUN = function(y) getM(object = object, x = y, untransform = space == "original", level = level))))
  if (Xclass == "factor")
   {
   	res[, 1] <- as.factor(res[, 1])
   	levels(res[, 1]) <- levels(mf[, 2])
   	rownames(res) <- NULL
   }
  if (Xclass != "factor")
   {
   	ord <- order(res[, 1])
   	res <- res[ord, ]
   	rownames(res) <- NULL
   }
  ymeasure <- switch(family(mod)$family, gaussian = "mean", binomial = "probability", poisson = "mean")
  Mname <- switch(family(mod)$family, gaussian = "mean(Y)", binomial = "P(Y)", poisson = "mean(Y)")
  Xname <- "X"
  if (space == "original" && family(mod)$family == "gaussian" && ypow != 1)
   {
    if (ypow == 0)
     {
      ymeasure <- "geometric mean"
      Mname <- "geomMean(Y)" 
      } else {
      ymeasure <- "median"
      Mname <- "median(Y)"
     }
    }
  if (space == "transformed")
   {
   	if (xpow != 1)
   	 {
   	  if (xpow == 0) Xname <- "log(X)" else Xname <- "X transf."
     }
    ymeasure <- switch(family(mod)$family, gaussian = "mean", binomial = "log(odds)", poisson = "log(mean)")
    Mname <- switch(family(mod)$family, gaussian = "mean(Y)", binomial = "log(odds(Y))", poisson = "log(mean(Y))")
    if (ypow != 1)
     {
      if (ypow == 0) Mname <- "mean(log(Y))" else Mname <- "mean(Ytrans)"
     }
    }
  names(res)[2] <- Mname
  if (Xclass != "factor")
   names(res)[1] <- Xname
  M <- list(M = res, ymeasure = ymeasure, space = space, ypow = ypow, xpow = xpow)
  class(M) <- "MY"
  M
 }
