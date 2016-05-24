"T2prob" <-
function(T) {
    if(any(T < 1)) {
      warning("Invalid return period")
      return(NULL)
    }
    return(1 - 1/T)
}
