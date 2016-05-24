# Get user parameter or default value from control object.
#
# @param control [\code{list}]\cr
#   Named list of control parameters.
# @param what [\code{character(1)}]\cr
#   Name of the parameter to get.
# @param default [\code{any}]\cr
#   Default value. This one is returned if there is not list element \dQuote{what}
#   in \dQuote{control}.
# @return [\code{any}]
getCMAESParameter = function(control, what, default) {
  return(coalesce(control[[what]], default))
}

norm2 = function(x) {
  return(drop(sqrt(crossprod(x))))
}
