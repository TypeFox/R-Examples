"cdfnor" <-
function(x,para) {
    if(! are.parnor.valid(para)) return()
    names(para$para) <- NULL # not needed but keep for auditing of package
    return(pnorm(x, mean = para$para[1], sd = para$para[2]))
}

