"quacau" <-
function(f,para,paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.parcau.valid(para)) return()
    }
    names(para$para) <- NULL
    return(qcauchy(f, location=para$para[1], scale=para$para[2]))
}
