se <-
function(x) {
  result <- sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))
  return(result)
}

