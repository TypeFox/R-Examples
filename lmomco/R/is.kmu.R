"is.kmu" <-
function(para) {
    if(para$type != "kmu") {
      warning("Parameters are not 'Kappa minus Mu' parameters")
      return(FALSE)
    }
    return(TRUE)
}

