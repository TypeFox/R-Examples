checkfornumber <-
function(x, name)
 {
  if ((!inherits(x, "integer") && !inherits(x, "numeric")) || length(x) > 1)
     stop(paste("if '", name, "' is provided, it must be a single number", sep = ""))
 }
