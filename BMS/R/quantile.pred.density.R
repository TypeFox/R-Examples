quantile.pred.density <-
function (x, probs = seq(0.25, 0.75, 0.25), names = TRUE, ...) 
{
    quout = .quantile.density(x$densities(), probs = probs, names = names, 
        normalize = FALSE)
    if (is.matrix(quout) && as.logical(names)) 
        rownames(quout) <- names(x$fit)
    return(quout)
}
