"is.aep4" <-
function(para) {
    if(para$type != "aep4") {
      warning("Parameters are not 4-p Asymmetric Exponential parameters.")
      return(FALSE)
    }
    return(TRUE)
}
