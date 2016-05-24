checkr <-
function(r)
 {
  if ((!inherits(r, "integer") && !inherits(r, "numeric")) || length(r) > 1)
       stop("the percent change in X, 'r', must be a number")
 }
