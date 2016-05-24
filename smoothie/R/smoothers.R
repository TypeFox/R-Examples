Fourier2d <-
function(x, bigdim=NULL, kdim=NULL) {
   ##
   ## Function to compute the Fast Fourier Transform (FFT) of a field 'x'.
   ##
   ## Arguments:
   ##
   ## 'x' 'n X m' numeric matrix.
   ## 'bigdim' numeric vector of length 2 giving the optimized dimension for the FFT.  If NULL, then 'kdim' must be supplied,
   ##   and it will be determined.
   ## 'kdim' dimension of the kernel matrix before expansion.  Only used if 'bigdim' is not supplied.
   ##
   ## Details: This function creates a matrix of zeros with twice the dimension of the input field, 'x'.
   ##   'x' is then put in the upper left corner of the matrix of zeros, missing values are set to zero,
   ##   and the FFT is calculated for the large matrix.  This is what is returned.  The output can then be
   ##   used by the 'kernel2dsmooth' function in order to possibly save an FFT in computing time at some point.
   ##
   ## See also: 'fft', 'kernel2dsmooth', 'hoods2d'
   ##
   ## Value: '2n X 2m' numeric matrix giving the Fourier transformed field.
   ##
   xdim <- dim( x)
   if( is.null( bigdim)) {
        if(is.null(kdim)) stop("Fourier2d: one of bigdim or kdim must be supplied.")
        bigdim <- xdim + kdim -1
        if( bigdim[1] <= 1024) bigdim[1] <- 2^ceiling(log2(bigdim[1]))
        else bigdim[1] <- ceiling( bigdim[1]/512)*512
        if( bigdim[2] <= 1024) bigdim[2] <- 2^ceiling(log2(bigdim[2]))
        else bigdim[2] <- ceiling( bigdim[2]/512)*512
   }
   out <- matrix( 0, bigdim[1], bigdim[2])
   out[1:xdim[1],1:xdim[2]] <- x
   out[ is.na( out)] <- 0
   return( fft( out))
} # end of 'Fourier2d' function.

kernel2dsmooth <-
function( x, kernel.type=NULL, K=NULL, W=NULL, X=NULL, xdim=NULL, Nxy=NULL, setup=FALSE, verbose=FALSE, ...) {
   ##
   ## Function to compute the boxcar kernel convolution smooth on a 2d image field.
   ##
   ## Arguments:
   ##
   ## 'x' 'k X r' matrix giving the grid point (pixel) values to be smoothed.
   ## 'kernel.type' (optional) character name of a kernel type allowed by 'kernel2dmeitsjer'.  One, and only one, of
   ##	this argument, 'K' or 'W' must be supplied.
   ## 'K' (optional) 'k1 X k2' matrix giving the kernel to be applied to the smoothing (e.g., if one wanted to
   ##   apply a different kernel other than the boxcar kernel, then these would be passed in here.
   ## 'W' (optional) '2k X 2r' matrix of Fourier transformed kernel weights.  If NULL, these will be found, and
   ##   if 'setup' is TRUE, these are what is returned, rather than a smoothed field.
   ## 'X' (optional) '2k X 2r' matrix giving the Fourier transformed expanded image of 'x'.  This is for when
   ##   multiple calls to this program are required (e.g., for different neighborhood lengths).  If supplied, it
   ##   will save one FFT in computation.  If 'W' is also supplied, then it will save two FFT's.
   ## 'xdim' (optional) numeric vector of length 2 giving the dimensions of 'x'.  If NULL, these are found here.
   ## 'Nxy' (optional) numeric giving the total number of grid points in 'x'.  If NULL, it will be calculated here.
   ## 'setup' logical, should just the Fourier transformed smoothing weights be calculated and returned?  If TRUE,
   ##   then 'W' should be NULL, or you will not get what you think.
   ## 'verbose' logical, should progress information be printed to the screen?
   ## '...' optional arguments to the 'kernel2dmeitsjer' function, only used if 'kernel.type' is supplied.
   ##
   ## Details:
   ##
   ## This 2-d spatial kernel smoother applies the smoother of Roberts and Lean (2008) to a spatial field.  Specifically,
   ## If X is a matrix of grid points, then the returned field, denoted by Ebert (2008) as <X>_s, is a smoothed field such
   ## that the value at each grid point '<X>_s[i,j]' is given by:
   ##
   ##   <X>_s[i,j] = sum_k sum_l X[i + k - 1 - (n-1)/2, j + l - 1 - (n-1)/2]*K[i + k - 1 - (n-1)/2, j + l - 1 - (n-1)/2],
   ##
   ## where n is the neighborhood length, k,l = 1, ..., n, and K[i + k - 1 - (n-1)/2, j + l - 1 - (n-1)/2] (default) is
   ## constant, and equal to 1/n^2.
   ##
   ## In order to be fast, loops are avoided completely.  Instead, the convolution theorem is applied with a Fast Fourier
   ## Transform (FFT).  If the weights 'W' are supplied, then you will save one FFT in computation time.  See Romanyuk (2005)
   ## for more information on this approach.
   ##
   ## The convolution theorem says that the Fourier transform of a convolution between two functions f and g is equal to the
   ## product of the Fourier transformed functions.  That is, if F denotes the Fourier transform, and * the convolution
   ## operator, F( f*g ) = k F(f)F(g), where 'k' is a scaling factor.  The neighborhood smooth is given by a convolution
   ## between the field and a boxcar kernel (i.e., a square around a point with constant value 1/n^2).  Because of the FFT,
   ## this enables a fast way to compute this convoultion.
   ##
   ## In order to zero-pad the field, and perform a cyclic convolution, it is necessary to expand the field, 'x', and re-arrange
   ## the boxcar kernel (or else it will not be centered on the points).
   ##
   ## References:
   ##
   ## Ebert E. E., 2008. Fuzzy verification of high resolution gridded forecasts: A review and proposed framework. Meteorol. Appl.,
   ##   15:51--64. DOI: 10.1002/met.25
   ##   Available at http://www.ecmwf.int/newsevents/meetings/workshops/2007/jwgv/METspecialissueemail.pdf
   ##
   ## Roberts, N. M. and H. W. Lean, 2008: Scale-selective verification of rainfall accumulations
   ##   from high-resolution forecasts of convective events.  Mon. Wea. Rev., 136:78--97.
   ##   DOI: 10.1175/2007MWR2123.1.
   ##
   ## Romanyuk, Y. A., 2005: "Discrete Signal Transformations."  Basics of Digital Signal Processing.  Part I.  Chapter 3.
   ##
   ## Value: If 'setup' is FALSE, then the smoothed field is returned, which is a numeric matrix with the same dimension as 'x'.
   ##   Otherwise, a '2k X 2r' matrix of Fourier transformed kernel weights is returned for future use with this function.
   ##
   # end of '!is.null(n)' stmt.
   if( is.null( xdim)) xdim <- dim( x)
   if( is.null( Nxy)) Nxy <- prod( xdim)
   if( is.null(W)) {
	if( !is.null(kernel.type)) K <- kernel2dmeitsjer( type=kernel.type, ...)
        else if( is.null(K)) stop("kernel2dsmooth: must give a value for at least one of kernel.type, K, or W")
	kdim <- dim( K)
        bigdim <- xdim + kdim - 1
        if( bigdim[1] <= 1024) bigdim[1] <- 2^ceiling(log2(bigdim[1]))
        else bigdim[1] <- ceiling( bigdim[1]/512)*512
        if( bigdim[2] <= 1024) bigdim[2] <- 2^ceiling(log2(bigdim[2]))
        else bigdim[2] <- ceiling( bigdim[2]/512)*512
        Kbig <- matrix( 0, bigdim[1], bigdim[2])
        kcen <- floor( (kdim+1)/2)
	Kbig[1:(kdim[1]-kcen[1]+1), 1:(kdim[2]-kcen[2]+1)] <- K[kcen[1]:kdim[1], kcen[2]:kdim[2]]
        if( kdim[1] > 1) Kbig[ (bigdim[1]-kcen[1]+2):bigdim[1], 1:(kdim[2]-kcen[2]+1)] <- K[1:(kcen[1]-1), kcen[2]:kdim[2]]
        if( kdim[2] > 1) Kbig[ 1:(kdim[1]-kcen[1]+1), (bigdim[2]-kcen[2]+2):bigdim[2]] <- K[ kcen[1]:kdim[1], 1:(kcen[2]-1)]
        if( all( kdim > 1)) Kbig[ (bigdim[1]-kcen[1]+2):bigdim[1], (bigdim[2]-kcen[2]+2):bigdim[2]] <- K[1:(kcen[1]-1), 1:(kcen[2]-1)]
        if( verbose) cat("Finding the FFT of the kernel matrix.\n")
        W <- fft( Kbig)/prod(bigdim)
        if( verbose) cat("FFT of kernel matrix found.\n")
        if( setup) return( W)
   } else bigdim <- dim( W)
   out <- matrix( 0, bigdim[1], bigdim[2])
   out[ 1:xdim[1], 1:xdim[2]] <- x
   out[ is.na( out)] <- 0
   # end of if 'is.null(W)' stmts.
   if( verbose) cat("Performing the convolution.\n")
   if( !is.null( X)) out <- Re( fft( X*W, inverse=TRUE))[1:xdim[1], 1:xdim[2]]
   else out <- Re( fft( fft( out)*W, inverse=TRUE))[1:xdim[1], 1:xdim[2]]
   if( verbose) cat("The convolution has been carried out.\n")
   return( zapsmall(out))
} # end of 'kernel2dsmooth' function.

kernel2dmeitsjer <- function(type="gauss", ...) {

   theta <- list(...)

   if(!is.null(theta$h)) h <- theta$h
   else if(!is.null(theta$nx) & !is.null(theta$ny)) {
        nx <- theta$nx
        ny <- theta$ny
        if(!is.null(theta$a)) a <- theta$a
        else a <- 1
        xgrid <- matrix( rep(1:nx,each=ny), nx, ny, byrow=TRUE)
        ygrid <- matrix( rep(1:ny,each=nx), nx, ny)
        xcen <- nx/2
        ycen <- ny/2
        h <- (xgrid - xcen)^2 + (ygrid - ycen)^2
   }

   if( !is.null( theta$sigma)) {
	sigma <- theta$sigma
	sigma2 <- sigma^2
   }

   if( type=="average") out <- matrix( 1/(theta$nx*theta$ny), theta$nx, theta$ny)
   else if( type=="boxcar") out <- matrix( 1/(theta$n^2), theta$n, theta$n)
   else if( type=="cauchy") out <- 1/(1+h/sigma)
   else if( type=="disk") {
	r   <- theta$r
 	r2 <- r^2
	rint  <- ceiling(r-0.5)
	diskg <- -rint:rint
	n <- length( diskg)
	x <- matrix( rep(diskg, each=n), n, n)
	y <- t(x)
	xy <- array( c(x, y), dim=c(n, n, 2))
	xy.absmax <- apply( xy, 1:2, function(x) max( abs(x))) # abs. distance from center.
	xy.absmin <- apply( xy, 1:2, function(x) min( abs(x))) # shortest direction to x- or y- coordinate at center.
	tmp1 <- tmp2 <- hold <- matrix(NA, n, n)
	val <- (xy.absmax+0.5)^2 + (xy.absmin-0.5)^2 # forming circular square centered at (0.5, 0.5)).
	id1 <- r2 < val # perimeter of circle
	id2 <- r2 >= val # interior of circle
	tmp1[id1] <- (xy.absmin[id1]-0.5)
	tmp1[id2] <- sqrt(r^2 - (xy.absmax[ id2] + 0.5)^2)
	val <- (xy.absmax-0.5)^2 + (xy.absmin+0.5)^2 # deal with points in corners.
	id1 <- r2 > val
	id2 <- r2 <= val
	tmp2[ id1] <- xy.absmin[id1]+0.5
	tmp2[id2] <- sqrt(r2 - (xy.absmax[id2] - 0.5)^2)
	val1 <- (xy.absmax+0.5)^2 + (xy.absmin+0.5)^2
	hold.id1 <- (( r2 < val1) & (r2 > (xy.absmax-0.5)^2 + (xy.absmin-0.5)^2))
	hold.id2 <- ((xy.absmin==0) & (xy.absmax-0.5 < r) & (xy.absmax+0.5>=r))
	hold.id <- hold.id1 | hold.id2
	hold[hold.id] <- r2*(0.5*(asin(tmp2[hold.id]/r) - asin(tmp1[hold.id]/r)) + 0.25*(sin(2*asin(tmp2[hold.id]/r)) - sin(2*asin(tmp1[hold.id]/r)))) - 
			(xy.absmax[hold.id]-0.5)*(tmp2[hold.id]-tmp1[hold.id]) + (tmp1[hold.id]-xy.absmin[hold.id]+0.5)
	hold[ is.na( hold)] <- 0
	hold <- hold + (val1 < r2)
	hold[rint+1,rint+1] <- min(pi*r2,pi/2)
	rc2 <- rint-0.5
     if ((rint>0) & (r > rc2) & (r2 < rc2^2+0.25)) { 
	tmp1  <- sqrt(r2 - rc2^2)
	tmp1n <- tmp1/r
	hold0 <- 2*(r^2*(0.5*asin(tmp1n) + 0.25*sin(2*asin(tmp1n)))-tmp1*(rint-0.5))
	hold[2*rint+1,rint+1]	<- hold0
	hold[rint+1,2*rint+1]	<- hold0
	hold[rint+1,1]		<- hold0
	hold[1,rint+1]		<- hold0
	hold[2*rint,rint+1]	<- hold[2*rint,rint+1] - hold0
	hold[rint+1,2*rint]	<- hold[rint+1,2*rint] - hold0
	hold[rint+1,2]		<- hold[rint+1,2]      - hold0
	hold[2,rint+1]		<- hold[2,rint+1]      - hold0
     } # end of if stmts.
     hold[rint+1,rint+1] <- min(hold[rint+1,rint+1],1)
     Sn <- sum( colSums( hold, na.rm=TRUE), na.rm=TRUE)
     out <- hold/Sn
   } else if(type=="epanechnikov") {
	out <- matrix(0, nx, ny)
	out[ h <= 1] <- 3/4*(1-h[ h <= 1]/sigma2)
   } else if( type=="exponential") {
	if( !is.null(theta$a)) a <- theta$a
	else a <- 1
	h <- sqrt( h)
	out <- a*exp(-h/(2*sigma2))
   } else if( type=="gauss") out <- (1/(2*pi*sigma2))*exp(-h/(2*sigma2))
   else if( type=="laplacian" | type=="unsharp") {
      if( !is.null( theta$a)) a <- theta$alpha
      else a <- 0
      a <- max( 0, min(a,1))
      o1 <- a/(a+1)
      o2 <- (1-a)/(1+a)
      out <- rbind( c(o1, o2, o1), c( o2, -4/(a+1), o2), c(o1, o2, o1))
      if( type=="unsharp") out <- rbind( rep(0,3), c(0,1,0), rep(0,3)) - out
   } else if(type=="LoG") {
      out <- (h-2*sigma2)*exp(-h/(2*sigma2))/(sigma2^2)
      out <- out - mean(out)
   } else if(type=="minvar") {
        out <- matrix(0, nx, ny)
	out[ h <= 1] <- 3/8*(3 - 5*h[ h <= 1]/sigma2)
   } else if( type=="multiquad") {
	a2 <- (theta$a)^2
	if( is.null( theta$inverse)) inverse <- FALSE
	else inverse <- theta$inverse
        out <- sqrt((h+a2))
	if( inverse) out <- 1/out
   } else if( type=="power") {
	p <- theta$p
        h <- sqrt( h)
	if( is.null( theta$do.log)) do.log <- FALSE
	else do.log <- theta$do.log
	if( !do.log) out <- -h^p
	else out <- -log(h^p+1)
   } else if( type=="prewitt") {
      out <- rbind( rep(1,3), rep(0,3), rep(-1,3))
      if( !is.null( theta$transpose)) if( theta$transpose) out <- t( out)
   } else if( type=="radial") {
	a <- theta$a
	m <- theta$m
	d <- theta$d
     	if( d%%2==0) out <- a*h^(2*m-d)*log(h)
     	else out <- a*h^(2*m-d)
     	out[ is.na(out)] <- 0
   } else if( type=="ratquad") {
	a <- theta$a
	out <- 1 - h/(h+a)
   } else if( type=="sobel") {
      out <- rbind( c(1,2,1), rep(0,3), c(-1,-2,-1))
      if( !is.null( theta$transpose)) if( theta$transpose) out <- t( out)
   } else if( type=="student") {
	p <- theta$p
	h <- sqrt(h)
	out <- 1/(1+h^p)
   } else if( type=="wave") {
	phi <- theta$phi
        h <- sqrt( h)
	out <- (phi/h)*sin(h/phi)
	out[ is.na( out)] <- 0
   }
   return(out)
} # end of 'kernel2dmeitsjer' function.

hoods2dsmooth <- function( x, lambda, W=NULL, setup=FALSE, ...) {
   if( floor( lambda) != lambda) {
        warning("hoods2dsmooth: attempting to give an illegal value for the neighborhood length.  Flooring lambda.")
        lambda <- floor( lambda)
   }
   if( lambda %% 2 == 0) {
	warning("hoods2dsmooth: attempting to give an even neighborhood length, subtracting one from it.")
	lambda <- lambda - 1
   }
   if( lambda < 1) {
	warning("hoods2dsmooth: attempting to give an illegal value for the neighborhood length.  Setting to one, and returning x untouched.")
	lambda <- 1
   }
   if( lambda == 1) return( x)
   else return( kernel2dsmooth( x=x, kernel.type="boxcar", n=lambda, W=W, setup=setup, ...))
} # end of 'hoods2dsmooth' function.

gauss2dsmooth <- function( x, lambda, W=NULL, setup=FALSE, ...) {

    out <- kernel2dsmooth( x=x, kernel.type="gauss", sigma=lambda, W=W, setup=setup, ...)
    return(out)

} # end of 'gauss2dsmooth' function.

disk2dsmooth <- function( x, lambda, W=NULL, setup=FALSE, ...) {
   return(kernel2dsmooth( x=x, kernel.type="disk", r=lambda, W=W, setup=setup, ...))
}

identity2dsmooth <- function( x, lambda=0, W=NULL, setup=FALSE, ...) {
   return(x)
}
