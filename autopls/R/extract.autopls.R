## extract.autopls.R Functions to extract values from an autopls object

predicted <- function (object)
{
  ## valid object?
  if (class (object) != 'autopls') stop ('needs object of class autopls')

  ## Get nlv
  lv <- get.lv (object)

  ## Get Yval
  Yval <- object$validation$pred [,,lv]
  return (Yval)
}

fitted.autopls <- function (object, ...)
{
  lv <- get.lv (object)
  class (object) <- 'mvr'
  object$fitted.values [,,lv]
}

coef.autopls  <- function (object, intercept = FALSE, ...)
{
  lv <- get.lv (object)
  class (object) <- 'mvr'
  rc <- coef (object, ncomp = lv, intercept = intercept) 
  nam <- rownames (rc)
  out <- as.vector (unlist (rc))
  names (out) <- nam
  out
}

coef.slim <- function (object, intercept = FALSE, ...)
{
  rc <- object$coefficients
  if (intercept) rc <- c(object$slimobj$intercept, rc)
  rc
} 

residuals.autopls <- function (object, ...)
{
  lv <- get.lv (object)
  class (object) <- 'mvr'
  residuals (object) [,,lv]
}

scores.autopls <- function (object, ...)
{
  lv <- get.lv (object)
  class (object) <- 'mvr'
  scores (object) [,lv]
}

loadings.autopls <- function (object, ...)
{
  lv <- get.lv (object)
  class (object) <- 'mvr'
  loadings (object) [,lv]
}

get.lv <- function (object)
{
  if (class (object) == 'autopls') return (object$metapls$current.lv)
  if (class (object) == 'slim') return (object$current.lv)
}

get.iter <- function (object)
{
  if (class (object) == 'autopls') return (object$metapls$current.iter)
  if (class (object) == 'slim') return (object$current.iter)
}

slim <- function (object)
{
  nobj <- list (coefficients = coef (object),
                method       = object$method,
                scale        = object$scale,
                call         = object$metapls$call,    
                predictors   = object$predictors,
                metapls = list (current.iter  = get.iter (object),
                           current.lv    = get.lv (object),
                           preprocessing = object$metapls$preprocessing,
                           scaling       = object$metapls$scaling,
                           val           = object$metapls$val),
                slimobj = list (intercept = coef (object, intercept = TRUE) [1], 
                           r2            = R2 (object, 'all'),
                           rmse          = RMSEP (object, 'all'),
                           N             = length (object$metapls$Y)))  
 
  class (nobj) <- 'slim'
  return (nobj)
}

