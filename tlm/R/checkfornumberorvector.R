checkfornumberorvector <-
function(x, name)
 {
  if ((!inherits(x, "integer") && !inherits(x, "numeric")))
     stop(paste("'", name, "' must be a number or a vector", sep = ""))
  }
