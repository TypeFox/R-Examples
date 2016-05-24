"ebayesthresh.wavelet.wd" <-
function (x.wd, vscale = "independent", smooth.levels = Inf,
    prior = "laplace", a = 0.5, bayesfac = FALSE, threshrule = "median")
{
    nlevs <- x.wd$nlevels
    slevs <- min(nlevs - 1, smooth.levels)
    if (is.character(vscale)) {
        vs <- substring(vscale, 1, 1)
        if (vs == "i")
            vscale <- mad(wavethresh::accessD(x.wd, level = nlevs - 1))
        if (vs == "l")
            vscale <- NA
    }
    for (j in (nlevs - slevs):(nlevs - 1)) {
        x.wd <-
            wavethresh::putD(x.wd, level = j,
                             v = ebayesthresh(wavethresh::accessD(x.wd, level = j),
                             prior, a, bayesfac, vscale, FALSE, threshrule))
    }
    return(x.wd)
}
