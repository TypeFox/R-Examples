"ebayesthresh.wavelet" <-
function (xtr, vscale = "independent", smooth.levels = Inf, 
    prior = "laplace", a = 0.5, bayesfac = FALSE, threshrule = "median")
{
    xcl <- class(xtr)
	if (class(xcl) == "dwt " && length(xcl) > 1) {
xtr <- ebayesthresh.wavelet.splus(xtr, vscale, smooth.levels, prior, a, bayesfac, threshrule)
return(xtr)}
	if (xcl == "wd") {
xtr <- ebayesthresh.wavelet.wd(xtr, vscale, smooth.levels, prior, a, bayesfac, threshrule)
return(xtr)}
if (xcl == "dwt"||xcl=="modwt") {
xtr <- ebayesthresh.wavelet.dwt(xtr, vscale, smooth.levels, prior, a, bayesfac, threshrule)
return(xtr)}
print("Unknown wavelet transform type; no smoothing performed")
return(xtr)
}
