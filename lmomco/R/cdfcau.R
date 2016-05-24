"cdfcau" <-
function(x,para) {
    if(! are.parcau.valid(para)) return()
    names(para$para) <- NULL
    return(pcauchy(x, location=para$para[1], scale=para$para[2]))
}

