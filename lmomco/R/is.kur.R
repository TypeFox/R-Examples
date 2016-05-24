"is.kur" <-
function(para) {
    if(para$type != "kur") {
      warning("Parameters are not Kumaraswamy parameters")
      return(FALSE)
    }
    return(TRUE)
}

