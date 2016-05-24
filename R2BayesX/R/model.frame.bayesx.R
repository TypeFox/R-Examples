model.frame.bayesx <- function (formula, ...) 
{
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  mn <- c("formula", "bayesx.setup", "bayesx.prg", "bayesx.run", "fitted.values", "residuals", 
    "response", "effects", "smooth.hyp", "random.hyp", "fixed.effects", "variance", "model.fit", 
    "call")
  check <- all("bayesx.setup" %in% names(formula))  
  mcheck <- NULL
  if(!check && is.list(formula))
    for(k in 1L:length(formula))
      mcheck <- c(mcheck, all("bayesx.setup" %in% names(formula[[k]])))
  if((length(nargs) || is.null(formula$bayesx.setup$data)) && !any(mcheck)) {
    fcall <- formula$call
    if(is.null(fcall)) {
      dots$formula <- formula
      rval <- do.call("parse.bayesx.input", dots)
      rval$data <- bayesx.reorder(rval$order, rval$data)
    } else {
      fcall[[1L]] <- as.name("parse.bayesx.input")
      fcall[names(nargs)] <- nargs
      env <- environment(terms(formula))
      if(is.null(env)) 
        env <- parent.frame()
      rval <- eval(fcall, env)
      rval$data <- bayesx.reorder(formula, rval$data)
    }
    if(!is.null(rval$h.random))
      return(gad(rval, rval$data))
    else
      return(rval$data)
  } else {
    if(any(mcheck)) {
      rval <- list()
      for(k in 1L:length(formula))
        if(mcheck[k])
          rval[[k]] <- model.frame.bayesx(formula[[k]], ...)
      if(length(rval) < 2L)
        rval <- rval[[1L]]
      return(rval)
    }
    if(!is.null(formula$bayesx.setup$h.random)) {
      rval <- gad(formula$bayesx.setup, formula$bayesx.setup$data)
      return(rval)
    } else {
      data <- bayesx.reorder(formula, formula$bayesx.setup$data)
      return(data)
    }
  }
}


gad <- function(x, dat)
{
  if(!is.null(x$h.random)) {
    dat <- list(dat)
    for(k in 1L:length(x$h.random)) {
      dat <- c(dat, list(x$h.random[[k]]$data))
      dat <- c(dat, gad(x$h.random[[k]], dat))
    }
    return(dat)
  } else return(NULL)
}

