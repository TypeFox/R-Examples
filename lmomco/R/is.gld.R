"is.gld" <-
function(para) {
    if(para$type != "gld") {
      warning("Parameters are not Generalized Lambda parameters")
      return(FALSE)
    }
    return(TRUE)
}

