"is.emu" <-
function(para) {
    if(para$type != "emu") {
      warning("Parameters are not 'Eta minus Mu' parameters")
      return(FALSE)
    }
    return(TRUE)
}

