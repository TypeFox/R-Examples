getComment <- function(x) {
  return(attr(x,"comment"))
}

"getComment<-" <- function(x,value) {
  attr(x,"comment")<-value
  return(x)
}
