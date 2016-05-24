"is.ln3" <-
function(para) {
    if(para$type != "ln3") {
      warning("Parameters are not log-Normal3 parameters")
      return(FALSE)
    }
    return(TRUE)
}

