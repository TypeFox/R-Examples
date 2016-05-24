"qualap" <-
function(f,para,paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.parlap.valid(para)) return()
    }
    XI <- para$para[1]
    A  <- para$para[2]

    x <- vector(mode="numeric", length=length(f))
    x[f <= 0.5] <- XI + A * log(2*   f[f <= 0.5])
    x[f >  0.5] <- XI - A * log(2*(1-f[f >  0.5]))
    names(f) <- NULL
    return(x)
}

