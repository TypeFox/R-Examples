
#------------------------------
# Main functions (out.globals)
#------------------------------

#' Get a parameter value from \code{solarius} models.
#'
#' @name modelPar
#' @rdname modelPar
#'
#' @param mod
#'    An object of \code{solarPolygenic}, \code{solarMultipoint} or \code{solarAssoc} classes. 
#'    See \code{\link{solarPolygenicClass}}, \code{\link{solarMultipointClass}} and \code{\link{solarAssocClass}}.
#' @param par
#'    A character, the parameter name.
#' @param ...
#'    Additional arguments.
#' @return 
#'    A value of the given parameter.
#'
#' @export
modelPar <- function(mod, par, ...)
{
  switch(par,  
    "cores" = modelParCores(mod),
    "CPUtime" = modelParCPUtime(mod, ...),
    "NumBatches" = modelParNumBatches(mod),
    stop(paste0("switch for `par` value"))
  )
}

#' @rdname modelPar
#'
#' @param format
#'    A character, the format of the time value.
#'    The default value is \code{"sec"}. The second possible value is \code{"POSIX"}.
#'    This argument is only for \code{modelParCPUtime} function.
#'
#' @export
modelParCPUtime <- function(mod, format = "sec", ...) 
{ 
  switch(class(mod)[1],
    "solarAssoc" = {
      t <- mod$assoc$tprofile$cputime.sec
      switch(format,
        "sec" = t,
        "POSIX" = format(.POSIXct(t, tz = "GMT"), "%H:%M:%S"),
        stop("swith error for `format`"))
    },
    stop("swith error for class of `mod`"))
}

#' @rdname modelPar
#' @export
modelParCores <- function(mod) 
{ 
  switch(class(mod)[1],
    "solarAssoc" = mod$assoc$cores,
    stop("swith error for class of `mod`"))
}

#' @rdname modelPar
#' @export
modelParNumBatches <- function(mod) 
{ 
  switch(class(mod)[1],
    "solarAssoc" = length(mod$assoc$solar$cmd),
    stop("swith error for class of `mod`"))
}














