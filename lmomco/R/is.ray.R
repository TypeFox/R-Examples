"is.ray" <-
function(para) {
    if(para$type != "ray") {
      warning("Parameters are not Rayleigh parameters")
      return(FALSE)
    }
    return(TRUE)
}

