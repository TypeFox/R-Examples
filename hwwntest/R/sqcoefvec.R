sqcoefvec <-
function (m0, filter.number = 10, family = "DaubLeAsymm", resolution = 4096, 
    stop.on.error = FALSE, plot.it=FALSE) 
{
    ecode <- 0
    sp <- support(filter.number = filter.number, family = family)
    psi <- draw.default(filter.number = filter.number, family = family, 
        resolution = resolution, plot.it = FALSE, enhance = FALSE)
    phi <- draw.default(filter.number = filter.number, family = family, 
        resolution = resolution, plot.it = FALSE, enhance = FALSE, scaling.function = TRUE)
    psif <- function(x, psi) {
        approx(x = psi$x, y = psi$y, xout = x, yleft = 0, yright = 0)$y
    }
    phif <- function(x, phi) {
        approx(x = phi$x, y = phi$y, xout = x, yleft = 0, yright = 0)$y
    }
    bigf <- function(x, m0, ll, psi, phi) {
        ans <- ((psif(x, psi = psi))^2) * phif((2^m0) * x - ll, 
            phi = phi)
    }
    lowl <- (2^m0) * sp$psi.lh - sp$phi.rh
    upl <- (2^m0) * sp$psi.rh - sp$phi.lh
    llvec <- lowl:upl
    ecvec <- rep(0, length(llvec))
    for (i in 1:length(llvec)) {
	    if (plot.it==TRUE)	{
		xx <- seq(from=sp$psi.lh, to=sp$psi.rh, length=100)
		plot(xx, bigf(xx, m0=m0, ll=llvec[i], psi=psi, phi=phi), type="l")
		scan()
		}
        ians <- integrate(bigf, lower = sp$psi.lh, upper = sp$psi.rh, 
            m0 = m0, ll = llvec[i], psi = psi, phi = phi, stop.on.error = stop.on.error)
        if (ians$message == "the integral is probably divergent") {
            cat(ians$message, "\n")
            ecvec[i] <- ians$value
        }
        else if (ians$message != "OK") {
            ecode <- 1
            return(list(ians = ians, ecode = ecode))
        }
        else ecvec[i] <- ians$value
    }
    ecvec <- 2^(m0/2) * ecvec
    l <- list(ll = llvec, ecoef = ecvec, m0 = m0, filter.number = filter.number, 
        family = family, ecode = ecode, ians = ians)
    l
}
