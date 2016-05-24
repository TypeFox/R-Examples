confint.bayesx <-
function(object, parm, level = 0.95, model = NULL, ...)
{
  object <- get.model(object, model)
  n <- length(object)
  if(n > 1L) {
    cf <- coef(object)
    medf <- edf(object)
    rval <- list()
    for(i in 1L:n) {
      if(!is.null(cf[[i]]))
        rval[[i]] <- confbayesx(cf[[i]], parm, level, medf[i], ...)
      else
        rval[[i]] <- NA
    }
    names(rval) <- names(object)
  } else {
    class(object) <- "bayesx"
    cf <- coef(object)
    if(!grepl("MCMC", object[[1L]]$model.fit$method))
      medf <- edf(object)
    else medf <- NA
    if(!is.null(cf))
      rval <- confbayesx(cf, parm, level, medf, ...)
    else
      rval <- NA
  }

  return(rval)
}

