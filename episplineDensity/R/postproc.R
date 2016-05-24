postproc <- function (out, postproc.controls = postproc.control ())
{
#
#           out: nloptr output
# postproc.opts: list of postprocessing options
#
if (any (names (postproc.controls) == "numevalpts"))
    numevalpts <- postproc.controls$numevalpts
else
    numevalpts <- 10000

par <- out$solution
#
# Now the postprocessing. Produce the estimate!
#
supportlow <- out$epiparameters$m0
supporthigh <- out$epiparameters$mN
#
# The "estpts" are the x's, a long sequence from low to high. Then
# compute the exponential epi-spline at each point.
#
estpts <- seq (supportlow, supporthigh, length=numevalpts)
expepi <- epispline (out$epiparameters, par, estpts, exponentiate=TRUE)
#
# If a "true" density is supplied, evaluate it at each of 
# estpts. It was a little different in the matlab code. Not implemented yet.
#
if (any (names (postproc.controls) == "trueDensity")) {
    d <- postproc.controls$trueDensity
    if (any (names (postproc.controls) == "trueDensityParams")) {
       params <- postproc.controls$trueDensityParams
       params$x <- out$x
       dens <- do.call (d, params)
    }
    else dens <- do.call (d, list (x = out$x))
#
# Compute MSE of true vs. epi
#
    out$mseepi <- trapz (estpts, (expepi - dens)^2)
}
#
# Integrate. This should be 1, presumably, but it might miss by a 
# little. Scale the y values so that the trapz() call produces exactly 1.
#

if(any (names (postproc.controls) == "normalize.to.1")) {
   out$orig.integral <- trapz (estpts, expepi)
   expepi <- expepi / out$orig.integral
   out$integral <- trapz (estpts, expepi)
}
#
# Pictures
#
if (any (names (postproc.controls) == "pic.types"))
    pic.types <- postproc.controls$pic.types
    if (!is.null (pic.types) && any (pic.types == "1")) {
        plot (estpts, expepi, type = "l", xlab = "X", ylab = "Density Estimate", main = "Epispline Density Estimate")
        points (out$x, rep (0, length (out$x)), col="red", pch=1, cex=1.5)
    }
#
# This is the real density estimate. Save it.
#
out$x.pts <- estpts
out$y.est <- expepi

return (out)
}
