myLeg <-
function (vec, arrond) {
  x <- vec
  lx <- length(x)
  if (lx < 3) 
    stop("not enought classes")
  res <- character(lx )
  res
  for (i in 1:(lx ))
  {res[i] <- paste(round(x[i],arrond),sep="")
  }
  res
}
