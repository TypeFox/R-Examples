"ebayesthresh.wavelet.wd" <-
function (x.wd, vscale = "independent", smooth.levels = Inf, 
    prior = "laplace", a = 0.5, bayesfac = FALSE, threshrule = "median") 
{
    nlevs <- x.wd$nlevels
    slevs <- min(nlevs - 1, smooth.levels)
	vv<-matrix(0,1,slevs)
    if (is.character(vscale)) {
        vs <- substring(vscale, 1, 1)

#for independent vscale, there is often very sparse fine level coefficients
#which estimates the noise variance as zero.  This attempts to remedy it.
        if (vs == "i"){
	for (i in (nlevs-slevs):(nlevs-1)){
            vscale <- mad(accessD(x.wd, level = i))
        if (vscale == 0) {
            vscale <- mad(accessD(x.wd, level = i), center = 0)
        }
	vv[nlevs-i]<-vscale
	}
	vscale<-vv[min(which(vv>0))]
}
        if (vs == "l") 
            vscale <- NA
    }
#print(vv)
    for (j in (nlevs - slevs):(nlevs - 1)) {
        x.wd <- putD(x.wd, level = j, v = ebayesthresh(accessD(x.wd, 
            level = j), prior, a, bayesfac, vscale, FALSE, threshrule))
    }
    return(x.wd)
}

