"ebayesthresh.wavelet.splus" <-
function (x.dwt, vscale = "independent", smooth.levels = Inf, 
    prior = "laplace", a = 0.5, bayesfac = FALSE, threshrule = "median") 
{
        nlevs <- attributes(x.dwt)$n.levels
        slevs <- min(nlevs, smooth.levels)
        if (is.character(vscale)) {
            vs <- substring(vscale, 1, 1)
            if (vs == "i") 
                vscale <- mad(x.dwt[[nlevs + 1]])
            if (vs == "l") 
                vscale <- NA
        }
        for (j in ((nlevs - slevs + 2):(nlevs + 1))) x.dwt[[j]] <- ebayesthresh(as.vector(x.dwt[[j]]), 
            prior, a, bayesfac, vscale, FALSE, threshrule)
        return(x.dwt)
    }
