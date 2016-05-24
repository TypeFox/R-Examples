checkq <-
function(q)
 {
  if ((!inherits(q, "integer") && !inherits(q, "numeric")) || length(q) > 1 || q <= 0)
       stop("the multiplicative change in X, 'q', must be a positive number")
 }
