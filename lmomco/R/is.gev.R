"is.gev" <-
function(para) {
    if(para$type != "gev") {
      warning("Parameters are not Generalized Extreme Value parameters")
      return(FALSE)
    }
    return(TRUE)
}

