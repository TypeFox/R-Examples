"is.cau" <-
function(para) {
    if(para$type != "cau") {
      warning("Parameters are not Cauchy parameters")
      return(FALSE)
    }
    return(TRUE)
}

