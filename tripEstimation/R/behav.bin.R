`behav.bin` <-
function (z, pimgs, weights = NULL)
{
    if (is.null(weights) && attr(pimgs, "Z"))
        weights <- c(diff(unclass(attr(pimgs, "Xtimes"))/3600))
    if (is.null(weights)) weights <- rep(1, length(pimgs))
    
    
    if (!(length(weights) == length(pimgs))) stop("length of weights do not match length of p-img list")
    if (nrow(z) != length(pimgs))
        stop("dimensions of chain do not match length of p-img list")
#dm <- dim(z)

 
    for (k in 1:length(weights)) {
        pimgs[[k]] <- bin.pimg(pimgs[[k]], t(z[k, 1:2, ]), w = weights[k])
    }
    pimgs
}

