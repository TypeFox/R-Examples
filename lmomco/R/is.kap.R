"is.kap" <-
function(para) {
    if(para$type != "kap") {
      warning("Parameters are not Kappa parameters.")
      return(FALSE)
    }
    return(TRUE)
}

