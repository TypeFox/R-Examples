dae <- function (y, times, parms, dy, res = NULL, func = NULL,
    method = c("mebdfi", "daspk", "radau", "gamd", "bimd"), ...)
{
    if (is.null(method))
        method <- "mebdfi"     
    else if (is.function(method) & !is.null(res))
	  if("res" %in% names(formals(method)))
        out <- method(y=y, times=times, parms=parms, dy=dy, res=res, ...)  
      else
		  stop("the solver cannot handle implicit equations")  
    else if (is.function(method) & !is.null(func))
        out <- method(y=y, times=times, parms=parms, func=func, ...) 
    else if(!is.null(res)) out <- switch(match.arg(method),
      mebdfi = mebdfi(y=y, times=times, parms=parms, dy=dy, res=res, ...),
      daspk  = daspk (y=y, times=times, parms=parms, dy=dy, res=res, ...),
      gamd   = stop("gamd cannot handle implicit equations"),
      bimd   = stop("bimd cannot handle implicit equations"),
      radau   = stop("radau5 cannot handle implicit equations"))
    else out <- switch(match.arg(method),
      mebdfi = mebdfi(y=y, times=times, parms=parms, dy=dy, func=func, ...),
      daspk  = daspk (y=y, times=times, parms=parms, dy=dy, func=func, ...),
      radau  = radau (y=y, times=times, parms=parms, func=func, ...),
      gamd   = gamd  (y=y, times=times, parms=parms, func=func, ...),
      bimd   = bimd  (y=y, times=times, parms=parms, func=func, ...))
    return(out)
}
