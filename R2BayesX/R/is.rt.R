is.rt <-
function(x)
{
  if(length(x)) {
    n <- length(x)
    isrt <- rep(FALSE, n)
    for(i in 1L:length(x)) {
      split <- splitme(x[i])
      if(length(split) > 2L) {
        if(resplit(split[1L:2L]) == c("r("))
          isrt[i] <- TRUE
      }
      if(grepl("\"rsps\"", x[i], fixed = TRUE))
        isrt[i] <- TRUE
      if(grepl("\"re\"", x[i], fixed = TRUE))
        isrt[i] <- TRUE
      if(grepl("\"ra\"", x[i], fixed = TRUE))
        isrt[i] <- TRUE
      if(grepl("\"random\"", x[i], fixed = TRUE))
        isrt[i] <- TRUE
    }
  } else isrt <- FALSE

  return(isrt)
}

