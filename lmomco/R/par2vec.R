"par2vec" <-
function(para,...) {
    if(length(para) != 0) {
      return(para$para)
    } else {
      warning("Empty para component, not lmomco para object?")
      return()
    }
}
