"is.exp" <-
function(para) {
    if(para$type != "exp") {
      warning("Parameters are not exponential parameters")
      return(FALSE)
    }
    return(TRUE)
}

