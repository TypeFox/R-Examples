"is.rice" <-
function(para) {
    if(para$type != "rice") {
      warning("Parameters are not Rice parameters")
      return(FALSE)
    }
    return(TRUE)
}

