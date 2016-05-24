# cartesian_product()
# Alternative to expand.grid() that takes the names of variables
# and environment in which to evaluate them.

cartesian_product <- function(variable.names, envir = parent.frame())
{
  size <- prod(vapply(variable.names,
                      function (variable) {length(get(variable, envir))},
					  1))
  res <- matrix(NA, ncol = length(variable.names), nrow = size)
  for (i in rev(1:length(variable.names)))
  {
    successor.size <- prod(sapply(variable.names[i:length(variable.names)],
                                  function (variable)
                                  {
                                    length(get(variable, envir))
                                  }))
    current.length <- length(get(variable.names[i], envir))
	# This fails for factors, which are not pure vectors.
	# Let's see if this solution works.
    res[, i] <- rep(as.vector(get(variable.names[i], envir)),
                    size / successor.size,
                    each = successor.size / current.length)
  }
  res <- data.frame(res, stringsAsFactors = FALSE)
  names(res) <- variable.names
  return(res)
}
