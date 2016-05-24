"pdfnor" <-
function(x,para) {
    if(! are.parnor.valid(para)) return()
    names(para$para) <- NULL
    return(dnorm(x, mean = para$para[1], sd = para$para[2]))
}

