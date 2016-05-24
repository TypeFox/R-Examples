checkInput <-
function(x, name, class, length, value, range, min, max, minEq, maxEq){

  ## check class (note: 'integer' is treated as dominant class)
  if (!missing(class))
    if (any(class == "integer")){
      if (any(class(x) != "integer") &
          (any(class(x) != "numeric") || any(x%%1 != 0)))
        stop(paste("'", name, "' must be of class integer", sep = ""))
    } else {
      if (!any(class(x) == class))
        stop(paste("'", name, "' must be of class ",
          paste(class, collapse = " OR "), sep = ""))
    }

  ## check length
  if (!missing(length))
    if (length(x) != length)
      stop(paste("'", name, "' must be of length ", length, sep = ""))

  ## check value
  if (!missing(value)){
    test <- logical()
    for (i in seq_along(x)) test[i] <- !any(value == x[i])
    if (any(test))
      stop(paste("'", name, "' cannot take values other than ",
                 paste(value, collapse = " OR "), sep = ""))
  }

  ## check range
  if (!missing(range))
    if (any(x < range[1]) || any(x > range[2]))
      stop(paste("'", name, "' cannot take values outside (",
                 range[1], ",", range[2], ")", sep = ""))

  ## check min
  if (!missing(min))
    if (any(x < min))
      stop(paste("'", name, "' cannot be smaller than ", min, sep = ""))

  ## check max
  if (!missing(max))
    if (any(x > max))
      stop(paste("'", name, "' cannot be larger than ", max, sep = ""))

  ## check maxEq
  if (!missing(maxEq))
    if (any(x >= maxEq))
      stop(paste("'", name, "' must be smaller than ", maxEq, sep = ""))

  ## check minEq
  if (!missing(minEq))
    if (any(x <= minEq))
      stop(paste("'", name, "' must be larger than ", minEq, sep = ""))
}
