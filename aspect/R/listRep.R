`listRep` <-
function(x,m) {
  y <- list()
  for (j in 1:m) y <- c(y,list(x))
  return(y)
}

