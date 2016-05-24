"quanor" <-
function(f,para,paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.parnor.valid(para)) return()
    }
    names(para$para) <- NULL
    return(qnorm(f,mean = para$para[1], sd = para$para[2]))
}

