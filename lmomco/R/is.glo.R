"is.glo" <-
function(para) {
    if(para$type != "glo") {
      warning("Parameters are not Generalized Logistic parameters")
      return(FALSE)
    }
    return(TRUE)
}

