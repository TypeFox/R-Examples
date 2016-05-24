"is.texp" <-
function(para) {
    if(para$type != "texp") {
      warning("Parameters are not truncated exponential parameters")
      return(FALSE)
    }
    return(TRUE)
}

