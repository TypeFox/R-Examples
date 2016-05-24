"is.gam" <-
function(para) {
    if(para$type != "gam") {
      warning("Parameters are not gamma parameters")
      return(FALSE)
    }
    return(TRUE)
}

