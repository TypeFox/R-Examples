"is.nor" <-
function(para) {
    if(para$type != "nor") {
      warning("Parameters are not Normal parameters")
      return(FALSE)
    }
    return(TRUE)
}

