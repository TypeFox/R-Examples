### Armin (first: 2012-09-21; Last: 2012-09-25)

### get the saturated model (with respect to an LV in focus)
saturate <- function(model, LV, ...)
  UseMethod("saturate", model)

saturate.plsm <- function(model, LV, ...){
  LVs <- model$latent
  Dn <- semPLS:::reorder(model$D)$Dn[, LV]
  allpred <- names(Dn)[Dn > 1]
  model <- semPLS:::addPath.plsm(model, from = allpred, to = LV)
  return(model)
}

### fit the saturated model (with respect to an LV in focus)
fitsaturated <- function(object, LV, ...)
  UseMethod("fitsaturated", object)

fitsaturated.sempls <- function(object, LV, ...){
  model_sat <- saturate(model = object$model, LV = LV)
  object$model <- model_sat
  objectFull <- semPLS:::resempls(object, data = object$data,
    method = "Standard")
  return(objectFull)
}
