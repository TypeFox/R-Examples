"is.st3" <-
function(para) {
    if(para$type != "st3") {
      warning("Parameters are not 3-parameter Student t parameters")
      return(FALSE)
    }
    return(TRUE)
}

