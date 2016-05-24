"is.tri" <-
function(para) {
    if(para$type != "tri") {
      warning("Parameters are not Triangular parameters")
      return(FALSE)
    }
    return(TRUE)
}

