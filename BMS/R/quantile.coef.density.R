quantile.coef.density <-
function (x, probs = seq(0.25, 0.75, 0.25), names = TRUE, ...) 
{
    quout = .quantile.density(x, probs = probs, names = names, 
        normalize = TRUE)
    if (is.matrix(quout) && as.logical(names)) 
        rownames(quout) <- sapply(x, function(lx) lx[["data.name"]])
    return(quout)
}
