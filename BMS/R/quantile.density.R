quantile.density <-
function (x, probs = seq(0.25, 0.75, 0.25), names = TRUE, normalize = TRUE, 
    ...) 
{
    my.quantile.density = function(x, probs, names, normalize, 
        ...) {
        ycs = (cumsum(x$y) - (x$y - x$y[[1]])/2) * diff(x$x[1:2])
        if (normalize) 
            ycs = ycs/(ycs[[length(ycs)]])
        xin = x$x
        maxi = length(ycs)
        qqs = sapply(as.list(probs), function(qu) {
            iii = sum(ycs <= qu)
            if (iii == maxi) 
                return(Inf)
            else if (iii == 0L) 
                return(-Inf)
            else {
                return(xin[[iii + 1]] + ((ycs[[iii + 1]] - qu)/(ycs[[iii + 
                  1]] - ycs[[iii]])) * (xin[[iii]] - xin[[iii + 
                  1]]))
            }
        })
        if (as.logical(names)) 
            names(qqs) = paste(format(100 * probs, trim = TRUE, 
                digits = max(2L, getOption("digits"))), "%", 
                sep = "")
        return(qqs)
    }
    probs = as.vector(probs)
    if (is.element("density", class(x))) 
        return(my.quantile.density(x = x, probs = probs, names = names, 
            normalize = normalize))
    if (!all(sapply(x, function(dd) is.element("density", class(dd))))) 
        stop("x needs to be a density or list of densities")
    if (length(x) == 1L) 
        return(my.quantile.density(x = x[[1]], probs = probs, 
            names = names, normalize = normalize))
    qout = sapply(x, my.quantile.density, probs = probs, names = FALSE, 
        normalize = normalize)
    if (!is.matrix(qout)) {
        if (length(probs) > 1) 
            return(qout)
        qout = as.matrix(qout)
    }
    else qout = t(qout)
    if (as.logical(names)) 
        colnames(qout) = paste(format(100 * probs, trim = TRUE, 
            digits = max(2L, getOption("digits"))), "%", sep = "")
    return(qout)
}
