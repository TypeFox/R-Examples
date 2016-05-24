cv <-
function(x,abs=TRUE,pc=TRUE) {
  result <- sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)
  if (abs==TRUE) {result <- abs(result)}
  if (pc==TRUE) {result <- 100*result}
  return(result)
}

