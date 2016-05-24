"is.gov" <-
function(para) {
    if(para$type != "gov") {
      warning("Parameters are not Govindarajulu parameters")
      return(FALSE)
    }
    return(TRUE)
}

