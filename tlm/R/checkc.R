checkc <-
function(c)
 {
  if ((!inherits(c, "integer") && !inherits(c, "numeric")) || length(c) > 1)
       stop("the additive change in X, 'c', must be a number")
 }
