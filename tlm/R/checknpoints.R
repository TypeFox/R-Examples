checknpoints <-
function(npoints)
 {
  if ((!inherits(npoints, "integer") && !inherits(npoints, "numeric")) || length(npoints) > 1 || ceiling(npoints) != floor(npoints) || npoints < 1)
       stop("'npoints' must be a positive integer")
 }
