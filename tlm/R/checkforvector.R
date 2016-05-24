checkforvector <-
function(x, name)
 {
  mess <- paste("'", name, "' must be a numeric vector", sep = "") 
  cond <- !inherits(x, "integer") && !inherits(x, "numeric")
  if (cond || (!cond && length(x) <= 1))
     stop(mess)
 }
