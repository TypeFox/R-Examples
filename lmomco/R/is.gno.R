"is.gno" <-
function(para) {
    if(para$type != "gno") {
      warning("Parameters are not Generalized Normal parameters")
      return(FALSE)
    }
    return(TRUE)
}

