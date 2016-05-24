effect <-
function(object, x1 = NULL, x2 = NULL, c = NULL, q = NULL, r = NULL, npoints = NULL, level = 0.95, nboot = 999, seed = 4321, verbose = TRUE)
 {  
  if (!inherits(object, "tlm"))
     stop("argument 'object' must be of class 'tlm'")
     
  if(any(is.na(coef(object$model))))
     stop("effect is not available for models with any missing estimated coefficient")
  
  if (!inherits(level, "numeric") || level <= 0 || level >= 1 || length(level) != 1)
    stop("'level' must be a number in (0, 1)")
  
  checknbootseed(nboot = nboot, seed = seed)
  if (verbose)
   cat("\nComputing effects...\n\n")
  set.seed(seed)
  modeltype <- modelType(object = object)     
  mod <- object$model
  mf <- model.frame(mod)
  mt <- attr(mf, "terms")
  dessignMatrix <- model.matrix(mt, data = mf)
  X <- colMeans(dessignMatrix)
  Xclass <- attr(mt, "dataClasses")[2] 
  if (Xclass == "factor")
   {
   	res <- effectFactor(object = object, level = level, nboot = nboot, modeltype = modeltype)
   	res$redux <- -1
   	unusedarg <- FALSE
   	if (!is.null(c(x1, x2, c, q, r, npoints)))
   	 unusedarg <- TRUE
   	res$unusedarg <- unusedarg
   	} else {
    res <- effectContinuous(object = object, x1 = x1, x2 = x2, c = c, q = q, r = r, npoints = npoints, level = level, nboot = nboot, modeltype = modeltype)
   }
  res$modeltype <- modeltype
  # check for the order of CI borders:
  effect <- res$effect
  CIlowername <- paste("lower", 100 * level, "%", sep = "")  
  # 1 or 2 effects?:
  neffects <- sum(colnames(effect) == CIlowername)
  aux <- effect
  for (i in 1:neffects)
   {
   	# where is the i-th "lower"?
   	wherelower <- which(colnames(effect) == CIlowername)[i]
   	if (any(effect[, wherelower] >= effect[, wherelower + 1]))
   	 {
   	  aux[, wherelower] <- effect[, wherelower + 1]
   	  aux[, wherelower + 1] <- effect[, wherelower]
   	 }
   }
  res$effect <- aux
  attr(res, "redux") <- res$redux
  res$redux <- NULL
  attr(res, "unusegarg") <- res$unusedarg
  res$unusegarg <- NULL
  attr(res, "modeltype") <- res$modeltype
  res$modeltype <- NULL  
  class(res) <- "effect"
  res
 }
