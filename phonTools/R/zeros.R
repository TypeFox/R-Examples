# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

zeros = function (x,y=0){
  if (missing(x)) return (0)
  if (length(x) == 1 & y == 0) return (rep (0, x))
  if (length(x) == 1 & y > 0)  return (matrix (0, x,y))

  nr = nrow(x)
  nc = ncol(x)

  if (is.null(nr) | is.null(nr)){
    output = rep(0, length(x))
    return (output)
  }
  output = matrix (0, nr ,nc)
  return (output)
}
