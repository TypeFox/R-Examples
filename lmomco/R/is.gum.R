"is.gum" <-
function(para) {
    if(para$type != "gum") {
      warning("Parameters are not Gumbel parameters")
      return(FALSE)
    }
    return(TRUE)
}

