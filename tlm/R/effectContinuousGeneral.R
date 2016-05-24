effectContinuousGeneral <-
function(object, x1, x2, c, q, r, npoints, level, nboot, modeltype)
 {
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
  # Default: x1 = Q1, x2 = Q3:
  IQR <- as.numeric(quantile(xobs, probs = c(1, 3) / 4))
  IQR <- data.frame(x1 = IQR[1], x2 = IQR[2])
  x1x2 <- NULL
  effectsOutRange <- F
  IQReffect <- F

  # Default case: x1 not provided:
  if (is.null(x1))
   {
   	IQReffect <- T
   	# Default:
   	x1x2 <- IQR
   }

  # If only x1 is provided:
  if (!is.null(x1) && is.null(x2) && is.null(c) && is.null(q) && is.null(r) && is.null(npoints))
   {
    cond <- !inherits(x1, "integer") && !inherits(x1, "numeric")
    if (cond || (!cond && length(x1) <= 1))
     stop("length of 'x1' must be > 1 if other arguments are not provided")
    l <- length(x1)
   	x1x2 <- data.frame(x1 = x1[-l], x2 = x1[-1])
   }

  # If x1 and x2 are provided:
  if (!is.null(x1) && !is.null(x2))
   {
    if (!inherits(x1, "integer") && !inherits(x1, "numeric"))
     stop("'x1' must be numeric")
    if (!inherits(x2, "integer") && !inherits(x2, "numeric"))
     stop("'x2' must be numeric")
    if (length(x1) != length(x2))
   	   stop("'x1' and 'x2' must have the same length")
    x1x2 <- data.frame(x1 = x1, x2 = x2)
   }

  # prevails c -> q -> r:

  # If providing x1 and c, but no x2, npoints:
  if (!is.null(x1) && is.null(x2) && !is.null(c) && is.null(npoints))
   {
    checkfornumberorvector(x = x1, name = "x1")
    checkc(c = c)
    x2 <- x1 + c
    x1x2 <- data.frame(x1 = x1, x2 = x2)
   }

  # If providing x1 number, c and npoints (no x2):
  if (!is.null(x1) && is.null(x2) && !is.null(c) && !is.null(npoints))
   {
    checkfornumber(x = x1, name = "x1")
    checkc(c = c)
    checknpoints(npoints = npoints)
    x1 <- x1 + 0:(npoints - 1) * c
    x2 <- x1 + c
    x1x2 <- data.frame(x1 = x1, x2 = x2)
   }

  # If providing x1 and q, but no x2, c, npoints:
  if (!is.null(x1) && is.null(x2) && is.null(c) && !is.null(q) && is.null(npoints))
   {
    checkfornumberorvector(x = x1, name = "x1")
    checkq(q = q)
    x2 <- x1 * q
    x1x2 <- data.frame(x1 = x1, x2 = x2)
   }

  # If providing x1 number, q and npoints, but no x2, c:
  if (!is.null(x1) && is.null(x2) && is.null(c) && !is.null(q) && !is.null(npoints))
   {
    checkfornumber(x = x1, name = "x1")
    checkq(q = q)
    checknpoints(npoints = npoints)
    x1 <- x1 * q ^ (0:(npoints - 1))
    x2 <- x1 * q
    x1x2 <- data.frame(x1 = x1, x2 = x2)
   }

  # If providing x1 and r, but no x2, c, q, npoints:
  if (!is.null(x1) && is.null(x2) && is.null(c) && is.null(q) && !is.null(r) && is.null(npoints))
   {
    checkfornumberorvector(x = x1, name = "x1")
    checkr(r = r)
    x2 <- x1 * (1 + r / 100)
    x1x2 <- data.frame(x1 = x1, x2 = x2)
   }

  # If providing x1 number, r and npoints, but no x2, c, q:
  if (!is.null(x1) && is.null(x2) && is.null(c) && is.null(q) && !is.null(r) && !is.null(npoints))
   {
    checkfornumber(x = x1, name = "x1")
    checkr(r = r)
    checknpoints(npoints = npoints)
    q <- 1 + r / 100
    x1 <- x1 * q ^ (0:(npoints - 1))
    x2 <- x1 * q
    x1x2 <- data.frame(x1 = x1, x2 = x2)
   }

  # Other possible cases:
  if (is.null(x1x2))
   {
   	IQReffect <- T
   	# Default:
   	x1x2 <- IQR
   }
   
  if (min(x1x2) < xmin || max(x1x2) > xmax)
   effectsOutRange <- T
  
  # No/Yes bootstrap cases (depending on the effect scale):
  if (any(modeltype %in% c(1, 3, 5, 6, 7, 9, 10, 12)))
   {
   	eval(parse(text = paste("res <- effectContinuousmod", modeltype, "(object = object, x1x2 = x1x2, level = level, nboot = nboot, IQReffect = IQReffect, effectsOutRange = effectsOutRange)", sep = "")))
   	} else {
    # Bootstrap cases:
    effect <- t(apply(x1x2, 1, FUN = function (x) geteffectx1x2(object = object, x1 = x[1], x2 = x[2], level = level, nboot = nboot)))
    effect <- cbind(x1x2, effect)
    names(effect) <- c(paste("x", 1:2, sep = ""), rep(c("Estimate", paste(c("lower", "upper"), 100 * level, "%", sep = "")), 2))
    names(effect)[c(3, 6)] <- paste(names(effect)[3], c("Diff", "Percent"), sep = "")
    if (IQReffect)
     {
   	  info <- paste("Adjusted change in the median of the response variable when the explanatory\nvariable changes from the 1st to the 3rd quartile (confidence intervals based\non ", nboot, " bootstrap\nsamples)", sep = "")
   	  } else {
      info <- paste("Adjusted change in the median of the response variable when the explanatory\nvariable changes from x1 to x2 (confidence intervals based on ", nboot, " bootstrap\nsamples)", sep = "")
    }
    if (any(MY(object)$M$Estimate < 0))
     info <- paste("WARNING: percent scale for effects not suitable because of negative values\nfor the adjusted median.\n\n", info, sep = "")
    if (effectsOutRange)
     info <- paste("WARNING: computing effects out of the observed range of X.\n\n", info, sep = "")          
    info <- paste(info, ":\n\n", sep = "")
    res <- list(effect = effect, info = info)
   }
  res	
 }
