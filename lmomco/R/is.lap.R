"is.lap" <-
function(para) {
    if(para$type != "lap") {
      warning("Parameters are not Laplace parameters")
      return(FALSE)
    }
    return(TRUE)
}

