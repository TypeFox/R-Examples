#simp.approx <- function(x) {
#  res <- .Fortran("dcalc",
#   as.double(x),
#   approx = double(1),
#   volume = double(1),        
#   numerator = double(1),
#   denominator = double(1),
#   gamma = double(1),
#   PACKAGE = "localdepth"
#   )
#  return(res)
#}
