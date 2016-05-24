"is.sla" <-
function(para) {
    if(para$type != "sla") {
      warning("Parameters are not Slash parameters")
      return(FALSE)
    }
    return(TRUE)
}

