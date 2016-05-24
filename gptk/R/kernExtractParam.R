kernExtractParam <-
function (kern, only.values=TRUE, untransformed.values=FALSE) {
  funcName <- paste(kern$type, "KernExtractParam", sep="")
  func <- get(funcName, mode="function")

  params <- func(kern, only.values=only.values, untransformed.values=untransformed.values)

  if ( any(is.nan(params)) )
    warning("Parameter has gone to NaN.")

  if ( "transforms" %in% names(kern) && (length(kern$transforms) > 0)
      && !untransformed.values )
    for ( i in seq(along=kern$transforms) ) {
      index <- kern$transforms[[i]]$index
      funcName <- optimiDefaultConstraint(kern$transforms[[i]]$type)
      func <- get(funcName$func, mode="function")
      if (funcName$hasArgs)
        params[index] <- func(params[index], "xtoa", kern$transformArgs[[i]])
      else
        params[index] <- func(params[index], "xtoa")
    }

  return (params)
}
