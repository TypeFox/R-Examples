"is.gpa" <-
function(para) {
    if(para$type != "gpa") {
      warning("Parameters are not Generalized Pareto parameters")
      return(FALSE)
    }
    return(TRUE)
}

