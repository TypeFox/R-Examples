"is.gep" <-
function(para) {
    if(para$type != "gep") {
      warning("Parameters are not Generalized Exponential Poisson parameters")
      return(FALSE)
    }
    return(TRUE)
}

