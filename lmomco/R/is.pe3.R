"is.pe3" <-
function(para) {
    if(para$type != "pe3") {
      warning("Parameters are not Pearson Type III parameters")
      return(FALSE)
    }
    return(TRUE)
}

