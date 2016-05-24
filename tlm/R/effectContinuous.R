effectContinuous <-
function(object, x1, x2, c, q, r, npoints, level, nboot, modeltype)
 {
  redux <- 0
  unusedarg <- FALSE  # if unusedarg == TRUE, the effect can't be calculated for the given arguments. Computing the default effect (only model is considered).
  # min and max X:
  mod <- object$model
  mf <- model.frame(mod)
  mt <- attr(mf, "terms")
  dessignMatrix <- model.matrix(mt, data = mf)
  xobs <- dessignMatrix[, 2]
  xpow <- object$xpow
  if (xpow != 1)
   xobs <- powerUntransform(xt = xobs, power = xpow)
  xmin <- min(xobs)
  xmax <- max(xobs)
  cmax <- xmax - xmin
  rmax <- 100 * (xmax / xmin - 1)
  # Exact effect (i.e., the model has an exact effect which depends only on both beta and the change x1 --> x2):
  # Case A: continuous X:
  #   Case A1: only the model is provided. Then we assume x1 --> x2 = IQR.
  # or
  #   Case A2: the model and the appropriate increase ('c', 'q' or 'r') are provided (and we ignore other extra provided arguments like number of points of a vector for x1 points)
  # Then, for cases A1 and A2, we report a single number (the exact effect) and its CI, for the given change x1 --> x2 (default is IQR).
  
  # Default: x1 = Q1, x2 = Q3:
  IQR <- as.numeric(quantile(xobs, probs = c(1, 3) / 4))
  cIQR <- IQR[2] - IQR[1]
  rIQR <- 100 * (IQR[2] / IQR[1] - 1)
  IQReffect <- F
  
  # provided: the model (1, 3, 7 or 10), and 'c':
  if (modeltype %in% c(1, 3, 7, 10) & !is.null(c))
   {
   	checkc(c = c)
   	effectsOutRange <- c > cmax
   	redux <- modeltype
   }
   
  # provided: only the model (1, 3, 7 or 10) (it does not matter if 'q' and/or 'r' are provided):
  # if (modeltype %in% c(1, 3, 7, 10) & is.null(c) & is.null(q) & is.null(r) & is.null(x1) & is.null(x2))
  if (modeltype %in% c(1, 3, 7, 10) & is.null(c) & is.null(x1) & is.null(x2))
   {
   	if (!is.null(c(q, r)))
   	 unusedarg <- TRUE
   	c <- cIQR
   	IQReffect <- T
   	effectsOutRange <- c > cmax
   	redux <- modeltype
   } 
     
  # provided: the model (5, 6, 9 or 12), and 'q':
  if (modeltype %in% c(5, 6, 9, 12) & !is.null(q))
   {
   	checkq(q = q)
   	r <- 100 * (q - 1)
   	effectsOutRange <- r > rmax
   	redux <- modeltype
   }

  # provided: the model (5, 6, 9 or 12), and 'r':
  if (modeltype %in% c(5, 6, 9, 12) & is.null(q) && !is.null(r))
   {
   	checkr(r = r)
    effectsOutRange <- r > rmax
  	redux <- modeltype
   }

  # provided: only the model (5, 6, 9 or 12) (it does not matter if 'c' is provided):
  if (modeltype %in% c(5, 6, 9, 12) & is.null(q) & is.null(r) & is.null(x1) & is.null(x2))
   {
   	if (!is.null(c))
   	 unusedarg <- TRUE
   	r <- rIQR
   	IQReffect <- T
   	effectsOutRange <- r > rmax
   	redux <- modeltype
   } 

  if (redux %in% c(1, 3, 7, 10))
   {
   	eval(parse(text = paste("res <- effectContinuousReduxmod", modeltype, "(object = object, c = c, IQReffect = IQReffect, effectsOutRange = effectsOutRange, level = level)", sep = "")))
   	} else {
     if (redux %in% c(5, 6, 9, 12))
      {
       eval(parse(text = paste("res <- effectContinuousReduxmod", modeltype, "(object = object, r = r, IQReffect  = IQReffect, effectsOutRange = effectsOutRange, level = level)", sep = "")))
       } else {
       res <- effectContinuousGeneral(object = object, x1 = x1, x2 = x2, c = c, q = q, r = r, npoints = npoints, level = level, nboot = nboot, modeltype = modeltype)
      }
    }
  res$redux <- redux
  res$unusedarg <- unusedarg
  res
 }
