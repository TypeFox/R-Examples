checknbootseed <-
function(nboot, seed)
 {
  if ((!inherits(nboot, "numeric") && !inherits(nboot, "integer")) || nboot <= 0 || length(nboot) != 1 || ceiling(nboot) != floor(nboot))
      stop("'nboot' must be a positive integer")
  
  if ((!inherits(seed, "numeric") && !inherits(seed, "integer")))
   	  stop("'seed' must be an integer")
 }
