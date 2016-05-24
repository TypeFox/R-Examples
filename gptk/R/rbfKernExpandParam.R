rbfKernExpandParam <-
function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  kern$inverseWidth <- params[1]	## linear domain params, i.e. untransformed inverse-width and signal variance
  kern$variance <- params[2]

  return (kern)
}
