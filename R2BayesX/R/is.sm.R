is.sm <-
function(x)
{
  if(length(x)) {
    n <- length(x)
    issm <- rep(FALSE, n)
    for(i in 1L:length(x)) {
      split <- splitme(x[i])
      if(length(split) > 2L) {
        if(resplit(split[1L:3L]) == "te(")
          issm[i] <- TRUE
        if(resplit(split[1L:2L]) %in% c("s(", "r(", "f("))
          issm[i] <- TRUE
      }
      if(length(split) > 3L)
        if(resplit(split[1L:3L]) == "sx(")
          issm[i] <- TRUE 
    }
  } else issm <- FALSE

  return(issm)
}

