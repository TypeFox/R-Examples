effNames <-
function (effects = NULL, loci = NULL, nloc = 1) 
{
    n <- paste(rep(noia::effectsNames[1], nloc), collapse = "")
    if (!is.null(effects)) {
        for (i in 1:length(loci)) {
            substr(n, loci[i], loci[i] + 1) <- effects[i]
        }
    }
    return(n)
}
