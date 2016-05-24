gaborPatch <- function(sf,
		theta = 0,
		rad = (theta * pi)/180,
		pc = 1,
		sigma = 1/6,
		psi = 0,
		grating = c("cosine", "sine"),
		npoints = 100,
		trim = 0, 
		trim.col = .5, ...)
{
    if(length(npoints) == 1) npoints <- rep(npoints, 2)
    if(length(sigma) == 1) sigma <- rep(sigma, 2)
    X <- ((1:npoints[1L])/ npoints[1L]) - .5
    Xm <- matrix(rep(X, npoints[2L]), npoints[2L])
    Ym <- t(Xm)
    Xt <- Xm * cos(rad) + Ym * sin(rad)
    Yt <- -Xm * sin(rad) + Ym * cos(rad)
    grating <- match.arg(grating)
    if(grating == "cosine"){
        wave <- pc * cos(2*pi*Xt*sf + psi)
    } else {
        wave <- pc * sin(2*pi*Xt*sf + psi)
    }
    gauss <- exp( -.5 *( (Xt/sigma[1L])^2 + (Yt/sigma[2L])^2 ))
    gabor <- wave * gauss
    gabor[gauss < trim] <- (trim.col*2 - 1)
    image(z = gabor, zlim = c(-1,1), col = gray((0:npoints[1L])/(npoints[1L])), 
    	axes = FALSE, asp = npoints[2L]/npoints[1L], ...)
    invisible(gabor)
}
