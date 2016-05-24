"ebayesthresh.wavelet.dwt" <-
function (x.dwt, vscale = "independent", smooth.levels = Inf, 
    prior = "laplace", a = 0.5, bayesfac = FALSE, threshrule = "median") 
{
        nlevs <- length(x.dwt)-1
        slevs <- min(nlevs, smooth.levels)
        if (is.character(vscale)) {
            vs <- substring(vscale, 1, 1)
            if (vs == "i") 
                vscale <- mad(x.dwt[[1]])
            if (vs == "l") 
                vscale <- NA
        }
        for (j in 1:slevs) {
            x.dwt[[j]] <- ebayesthresh(x.dwt[[j]], 
                prior, a, bayesfac, vscale, FALSE, 
                threshrule)
        }
     return(x.dwt)
}
