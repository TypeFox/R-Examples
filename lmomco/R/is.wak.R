"is.wak" <-
function(para) {
    if(para$type != "wak") {
      warning("Parameters are not Wakeby parameters")
      return(FALSE)
    }
    return(TRUE)
}

