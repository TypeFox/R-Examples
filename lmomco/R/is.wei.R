"is.wei" <-
function(para) {
    if(para$type != "wei") {
      warning("Parameters are not Weibull parameters")
      return(FALSE)
    }
    return(TRUE)
}

