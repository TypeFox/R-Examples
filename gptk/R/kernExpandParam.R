kernExpandParam <-
function (kern, params, untransformed.values=FALSE) {
# browser()
  if ( is.list(params) )
    params <- params$values
  
  if ( "transforms" %in% names(kern) && (length(kern$transforms) > 0)
      && !untransformed.values )
    for ( i in seq(along=kern$transforms) ) {
      index <- kern$transforms[[i]]$index
      funcName <- optimiDefaultConstraint(kern$transforms[[i]]$type)
      func <- get(funcName$func, mode="function")
      if (funcName$hasArgs)
        params[index] <- func(params[index], "atox", kern$transformArgs[[i]]) ## log-transformed params just been exp-transformed
      else {
        params[index] <- func(params[index], "atox")
      }
    }

  funcName <- paste(kern$type, "KernExpandParam", sep="")
  func <- get(funcName, mode="function")
# browser()
  kern <- func(kern, params)

  return (kern)
  
}
