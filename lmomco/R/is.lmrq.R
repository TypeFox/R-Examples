"is.lmrq" <-
function(para) {
    if(para$type != "lmrq") {
      warning("Parameters are not Linear Mean Residual Quantile parameters")
      return(FALSE)
    }
    return(TRUE)
}

