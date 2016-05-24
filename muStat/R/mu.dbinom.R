`mu.dbinom` <-
function(x, size, prob, log=FALSE)
  if (size==0) 1 else dbinom(x, size, prob, log=FALSE)
