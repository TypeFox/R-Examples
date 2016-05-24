"is.revgum" <-
function(para) {
    if(para$type != "revgum") {
      warning("Parameters are not Reverse Gumbel parameters")
      return(FALSE)
    }
    return(TRUE)
}

