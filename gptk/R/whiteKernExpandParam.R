whiteKernExpandParam <-
function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  kern$variance <- params[1]	## linear domain param, i.e. the untransformed noise variance

  return (kern)
}
