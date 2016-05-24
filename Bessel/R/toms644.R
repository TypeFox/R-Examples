#### At first, simple interfaces to the main subroutines
#### in ../src/zbsubs.f
####	~~~~~~~~~~~~~~~
##   10: zbesh(zr, zi, fnu, kode, m, n, cyr, cyi, nz, ierr)

##  360: zbesi(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
##  631: zbesj(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
##  899: zbesk(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
## 1182: zbesy(zr, zi, fnu, kode, n, cyr, cyi, nz, cwrkr, cwrki, ierr)
##
## 1467: zairy(zr, zi,	id, kode, air, aii, nz, ierr)
## 1862: zbiry(zr, zi,	id, kode, bir, bii, ierr)

##  360: zbesi(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
BesselI <- function(z, nu, expon.scaled = FALSE, nSeq = 1)
{
    nz <- length(z)
    if(nz == 0) return(z)
    isNum <- is.numeric(z)
    if(isNum) {
	zr <- as.double(z)
	zi <- numeric(nz)
    } else if(is.complex(z)) {
	zr <- Re(z)
	zi <- Im(z)
    } else stop("'z' must be complex or numeric")
    nu <- as.double(nu)
    stopifnot(length(nu) == 1,
	      nSeq >= 1, nSeq == as.integer(nSeq))

    if(nu < 0) {
	## I(-nu,z) = I(nu,z) + (2/pi)*sin(pi*nu)*K(nu,z)  [ = A.& S. 9.6.2, p.375 ]
	if(nu == round(nu)) ## <==> sin(pi*nu) == 0
	    return(BesselI(z, -nu, expon.scaled, nSeq=nSeq))
        ## else
	nu. <- -nu + seq_len(nSeq) - 1
	kf <- rep(2/pi*sin(pi*nu.)* (if(expon.scaled) exp(-2*abs(z)) else 1),
		  each=nz)
	return(	  BesselI(z, -nu, expon.scaled, nSeq=nSeq) +
	       kf*BesselK(z, -nu, expon.scaled, nSeq=nSeq))
    }
    ## else  nu >= 0 :

    r <- if(isNum) numeric(nz * nSeq) else complex(nz * nSeq)
    if(nSeq > 1) r <- matrix(r, nz, nSeq)
    for(i in seq_len(nz)) {
	ri <- .Fortran("zbesi",
		       zr[i], zi[i],
		       fnu = nu,
		       kode= as.integer(1L + as.logical(expon.scaled)),
					# 1 or 2, exactly as desired
		       n = as.integer(nSeq),
		       cyr = double(nSeq),
		       cyi = double(nSeq),
		       nz   = integer(1),
		       ierr = integer(1), PACKAGE = "Bessel")
	if(ri$ierr) {
	    f.x <- sprintf("'zbesi(%g %s %gi, nu=%g)'", zr[i],
                           c("-","+")[1+(z >= 0)], abs(zi[i]), nu)
	    if(ri$ierr == 3)
		warning(sprintf("%s' large arguments -> precision loss\n (of at least half machine accuracy)", f.x))
	    else if(ri$ierr == 2) {
		if(getOption("verbose"))
                    message(sprintf("%s  -> overflow ; returning Inf\n", f.x))
                ri$cyr <- ri$cyi <- Inf
            }
	    else stop(sprintf("%s [Fortran] error ierr = %d", f.x, ri$ierr))
	}
	rz <- if(isNum) ri$cyr else complex(real = ri$cyr,
					    imaginary = ri$cyi)
	if(nSeq > 1) r[i,] <- rz else r[i] <- rz
    }
    r
} ## I()

##  631: zbesj(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
BesselJ <- function(z, nu, expon.scaled = FALSE, nSeq = 1)
{
    nz <- length(z)
    if(nz == 0) return(z)
    isNum <- is.numeric(z)
    if(isNum) {
	zr <- as.double(z)
	zi <- numeric(nz)
    } else if(is.complex(z)) {
	zr <- Re(z)
	zi <- Im(z)
    } else stop("'z' must be complex or numeric")
    nu <- as.double(nu)
    stopifnot(length(nu) == 1,
	      nSeq >= 1, nSeq == as.integer(nSeq))
    if(nu < 0) {
	## J(-fnu,z) = J(fnu,z)*cos(pi*fnu) - Y(fnu,z)*sin(pi*fnu)
	if(expon.scaled)
	    stop("'expon.scaled=TRUE' not yet implemented for nu < 0")
	pnu <- rep(pi*(-nu + seq_len(nSeq) - 1), each=nz)
	return(BesselJ(z, -nu, nSeq=nSeq)*cos(pnu) -
	       if(nu == round(nu)) 0 else
	       BesselY(z, -nu, nSeq=nSeq)*sin(pnu))
    }
    ## else  nu >= 0 :

    r <- if(isNum) numeric(nz * nSeq) else complex(nz * nSeq)
    if(nSeq > 1) r <- matrix(r, nz, nSeq)
    for(i in seq_len(nz)) {
	ri <- .Fortran("zbesj",
		       zr[i], zi[i],
		       fnu = nu,
		       kode= as.integer(1L + as.logical(expon.scaled)),
					# 1 or 2, exactly as desired
		       n = as.integer(nSeq),
		       cyr = double(nSeq),
		       cyi = double(nSeq),
		       nz   = integer(1),
		       ierr = integer(1), PACKAGE = "Bessel")
	if(ri$ierr) {
	    f.x <- sprintf("'zbesj(%g %s %gi, nu=%g)'", zr[i],
                           c("-","+")[1+(z >= 0)], abs(zi[i]), nu)
	    if(ri$ierr == 3)
		warning(sprintf("%s large arguments -> precision loss\n	 (of at least half machine accuracy)", f.x))
	    else if(ri$ierr == 2) {
		if(getOption("verbose"))
                    message(sprintf("%s  -> overflow ; returning Inf\n", f.x))
                ri$cyr <- ri$cyi <- Inf
            }
	    else stop(sprintf("%s [Fortran] error ierr = %d", f.x, ri$ierr))
	}
	rz <- if(isNum) ri$cyr else complex(real = ri$cyr,
					    imaginary = ri$cyi)
	if(nSeq > 1) r[i,] <- rz else r[i] <- rz
    }
    r
} ## J()

##  899: zbesk(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
BesselK <- function(z, nu, expon.scaled = FALSE, nSeq = 1)
{
    nz <- length(z)
    if(nz == 0) return(z)
    isNum <- is.numeric(z)
    if(isNum) {
	zr <- as.double(z)
	zi <- numeric(nz)
    } else if(is.complex(z)) {
	zr <- Re(z)
	zi <- Im(z)
    } else stop("'z' must be complex or numeric")
    nu <- as.double(nu)
    stopifnot(length(nu) == 1,
	      nSeq >= 1, nSeq == as.integer(nSeq))

    if(nu < 0) {
	## K(-nu,z) = K(nu,z)
	return(BesselK(z, -nu, expon.scaled, nSeq=nSeq))
    }
    ## else  nu >= 0 :

    r <- if(isNum) numeric(nz * nSeq) else complex(nz * nSeq)
    if(nSeq > 1) r <- matrix(r, nz, nSeq)
    for(i in seq_len(nz)) {
	ri <- .Fortran("zbesk",
		       zr[i], zi[i],
		       fnu = nu,
		       kode= as.integer(1L + as.logical(expon.scaled)),
					# 1 or 2, exactly as desired
		       n = as.integer(nSeq),
		       cyr = double(nSeq),
		       cyi = double(nSeq),
		       nz   = integer(1),
		       ierr = integer(1), PACKAGE = "Bessel")
	if(ri$ierr) {
	    f.x <- sprintf("'zbesk(%g %s %gi, nu=%g)'", zr[i],
                           c("-","+")[1+(z >= 0)], abs(zi[i]), nu)
	    if(ri$ierr == 3)
		warning(sprintf("%s large arguments -> precision loss\n	 (of at least half machine accuracy)", f.x))
	    else if(ri$ierr == 2) {
		if(getOption("verbose"))
                    message(sprintf("%s  -> overflow ; returning Inf\n", f.x))
                ri$cyr <- ri$cyi <- Inf
            }
	    else stop(sprintf("%s [Fortran] error ierr = %d", f.x, ri$ierr))
	}
	rz <- if(isNum) ri$cyr else complex(real = ri$cyr,
					    imaginary = ri$cyi)
	if(nSeq > 1) r[i,] <- rz else r[i] <- rz
    }
    r
} ## K()


## 1182: zbesy(zr, zi, fnu, kode, n, cyr, cyi, nz, cwrkr, cwrki, ierr)
BesselY <- function(z, nu, expon.scaled = FALSE, nSeq = 1)
{
    nz <- length(z)
    if(nz == 0) return(z)
    isNum <- is.numeric(z)
    if(isNum) {
	zr <- as.double(z)
	zi <- numeric(nz)
    } else if(is.complex(z)) {
	zr <- Re(z)
	zi <- Im(z)
    } else stop("'z' must be complex or numeric")
    nu <- as.double(nu)
    stopifnot(length(nu) == 1,
	      nSeq >= 1, nSeq == as.integer(nSeq))

    if(nu < 0) {
	## Y(-fnu,z) = Y(fnu,z)*cos(pi*fnu) + J(fnu,z)*sin(pi*fnu)
	if(expon.scaled)
	    stop("'expon.scaled=TRUE' not yet implemented for nu < 0")
	pnu <- rep(pi*(-nu + seq_len(nSeq) - 1), each=nz)
	return(BesselY(z, -nu, nSeq=nSeq)*cos(pnu) +
	       if(nu == round(nu)) 0 else
	       BesselJ(z, -nu, nSeq=nSeq)*sin(pnu))
    }
    ## else  nu >= 0 :

    r <- if(isNum) numeric(nz * nSeq) else complex(nz * nSeq)
    if(nSeq > 1) r <- matrix(r, nz, nSeq)
    for(i in seq_len(nz)) {
        if(zr[i] == 0 && zi[i] == 0) {
            rz <- if(isNum) -Inf else
            ## A limit  z -> 0  depends on the *direction*; it is
            ## ``a version of Inf"', i.e. the only *complex* Inf, if you think
            ## of the complex sphere. --> we use the same as 1/(0+0i):
            1/(0+0i)
        } else {
## 1182: zbesy(zr, zi, fnu, kode, n, cyr, cyi, nz, cwrkr, cwrki, ierr)
            ri <- .Fortran("zbesy",
                           zr[i], zi[i],
                           fnu = nu,
                           kode= as.integer(1L + as.logical(expon.scaled)),
					# 1 or 2, exactly as desired
                           n = as.integer(nSeq),
                           cyr = double(nSeq),
                           cyi = double(nSeq),
                           nz   = integer(1),
                           cwrkr= double(nSeq),
                           cwrki= double(nSeq),
                           ierr = integer(1), PACKAGE = "Bessel")
            if(ri$ierr) {
                f.x <- sprintf("'zbesy(%g %s %gi, nu=%g)'", zr[i],
                               c("-","+")[1+(z >= 0)], abs(zi[i]), nu)
                if(ri$ierr == 3)
                    warning(sprintf("%s large arguments -> precision loss\n	 (of at least half machine accuracy)", f.x))
                else if(ri$ierr == 2) {
                    if(getOption("verbose"))
                        message(sprintf("%s  -> overflow ; returning Inf\n", f.x))
                    ri$cyr <- ri$cyi <- Inf
                }
                else stop(sprintf("%s [Fortran] error ierr = %d", f.x, ri$ierr))
            }
            rz <- if(isNum) ri$cyr else complex(real = ri$cyr,
                                                imaginary = ri$cyi)
        }
        if(nSeq > 1) r[i,] <- rz else r[i] <- rz
    }
    r
} ## Y()

##---------------- Hankel function H() ------------------

##   10: zbesh(zr, zi, fnu, kode, m, n, cyr, cyi, nz, ierr)
BesselH <- function(m, z, nu, expon.scaled = FALSE, nSeq = 1)
{
## c***keywords	 H-bessel functions,bessel functions of complex argument,
## c		 bessel functions of third kind,hankel functions
## c***author  amos, donald e., sandia national laboratories
## c***purpose	to compute the H-bessel functions of a complex argument
## c***description
## c
## c			  ***a double precision routine***
## c	     on kode=1, zbesh computes an n member sequence of complex
## c	     hankel (bessel) functions	cy(j) = H(m,fnu+j-1,z) for kinds m=1 or 2,
## c	  real, nonnegative orders fnu+j-1, j=1,...,n, and complex
## c	     z != cmplx(0.,0.) in the cut plane -pi < arg(z) <= pi.
## c	     on kode=2, zbesh returns the scaled Hankel functions
## c
## c	     cy(i)= exp(-mm*z*i) * H(m,fnu+j-1,z)	mm=3-2*m,   i**2=-1.
## c
## c	     which removes the exponential behavior in both the upper and
## c	     lower half planes. definitions and notation are found in the
## c	     nbs handbook of mathematical functions (ref. 1).

    m <- as.integer(m)
    stopifnot(length(m) == 1, m == 1 || m == 2)
    nz <- length(z)
    if(nz == 0) return(z)
    isNum <- is.numeric(z)
    if(isNum) {
	zr <- as.double(z)
	zi <- numeric(nz)
    } else if(is.complex(z)) {
	zr <- Re(z)
	zi <- Im(z)
    } else stop("'z' must be complex or numeric")
    nu <- as.double(nu)
    stopifnot(length(nu) == 1,
	      nSeq >= 1, nSeq == as.integer(nSeq))

    if(nu < 0) {
	## H(1,-fnu,z) = H(1,fnu,z)*cexp( pi*fnu*i)
	## H(2,-fnu,z) = H(2,fnu,z)*cexp(-pi*fnu*i) ; i^2=-1
	if(expon.scaled)
	    stop("'expon.scaled=TRUE' not yet implemented for nu < 0")
	pnu <- rep(c(1i,-1i)[m] * pi*(-nu + seq_len(nSeq) - 1), each=nz)
	return(BesselH(m, z, -nu, nSeq=nSeq)* exp(pnu))
    }
    ## else  nu >= 0 :

    r <- if(isNum) numeric(nz * nSeq) else complex(nz * nSeq)
    if(nSeq > 1) r <- matrix(r, nz, nSeq)
    for(i in seq_len(nz)) {
	## zbesh(zr, zi, fnu, kode, m, n, cyr, cyi, nz, ierr)
	ri <- .Fortran("zbesh",
		       zr[i], zi[i],
		       fnu = nu,
		       kode= as.integer(1L + as.logical(expon.scaled)),
					# 1 or 2, exactly as desired
		       m = m,
		       n = as.integer(nSeq),
		       cyr = double(nSeq),
		       cyi = double(nSeq),
		       nz   = integer(1),
		       ierr = integer(1), PACKAGE = "Bessel")
	if(ri$ierr) {
	    f.x <- sprintf("'zbesh(%g %s %gi, nu=%g)'", zr[i],
                           c("-","+")[1+(z >= 0)], abs(zi[i]), nu)
	    if(ri$ierr == 3)
		warning(sprintf("%s large arguments -> precision loss\n	 (of at least half machine accuracy)", f.x))
	    else if(ri$ierr == 2) {
		if(getOption("verbose"))
                    message(sprintf("%s  -> overflow ; returning Inf\n", f.x))
                ri$cyr <- ri$cyi <- Inf
            }
	    else stop(sprintf("%s [Fortran] error ierr = %d", f.x, ri$ierr))
	}
	rz <- if(isNum) ri$cyr else complex(real = ri$cyr,
					    imaginary = ri$cyi)
	if(nSeq > 1) r[i,] <- rz else r[i] <- rz
    }
    r
} ## H()

##---------------- Airy Functions Ai()	Bi() ------------------

## 1471:  zairy(zr, zi, id, kode, air, aii, nz, ierr)
## 1867:  zbiry(zr, zi, id, kode, bir, bii, ierr)

AiryA <- function(z, deriv = 0, expon.scaled = FALSE)
{
## Purpose:  to compute airy functions ai(z) and dai(z) for complex z
##	     airy function : "Bessel functions of order one third"
    deriv <- as.integer(deriv)
    stopifnot(length(deriv) == 1, deriv == 0 || deriv == 1)
    nz <- length(z)
    if(nz == 0) return(z)
    isNum <- is.numeric(z)
    if(isNum) {
	zr <- as.double(z)
	zi <- numeric(nz)
    } else if(is.complex(z)) {
	zr <- Re(z)
	zi <- Im(z)
    } else stop("'z' must be complex or numeric")
    r <- if(isNum) numeric(nz) else complex(nz)
    for(i in seq_len(nz)) {
	## zairy(zr, zi, id, kode, air, aii, nz, ierr)
	ri <- .Fortran("zairy",
		       zr[i], zi[i],
		       id = deriv,
		       kode= as.integer(1L + as.logical(expon.scaled)),
					# 1 or 2, exactly as desired
		       air = double(1),
		       aii = double(1),
		       nz   = integer(1),
		       ierr = integer(1), PACKAGE = "Bessel")
	if(ri$ierr) {
	    f.x <- sprintf("'zairy(%g %s %gi, deriv=%d)'", zr[i],
                           c("-","+")[1+(z >= 0)], abs(zi[i]), deriv)
	    if(ri$ierr == 3)
		warning(sprintf("%s large arguments -> precision loss\n	 (of at least half machine accuracy)", f.x))
	    else if(ri$ierr == 2) {
		if(getOption("verbose"))
                    message(sprintf("%s  -> overflow ; returning Inf\n", f.x))
                ri$air <- ri$aii <- Inf
            }
	    else stop(sprintf("%s [Fortran] error ierr = %d", f.x, ri$ierr))
	}
	r[i] <- if(isNum) ri$air else complex(real = ri$air, imaginary = ri$aii)
    }
    r
} ## AiryA()

AiryB <- function(z, deriv = 0, expon.scaled = FALSE)
{
## Purpose:  to compute airy functions bi(z) and dbi(z) for complex z
##	     airy function : "Bessel functions of order one third"
    deriv <- as.integer(deriv)
    stopifnot(length(deriv) == 1, deriv == 0 || deriv == 1)
    nz <- length(z)
    if(nz == 0) return(z)
    isNum <- is.numeric(z)
    if(isNum) {
	zr <- as.double(z)
	zi <- numeric(nz)
    } else if(is.complex(z)) {
	zr <- Re(z)
	zi <- Im(z)
    } else stop("'z' must be complex or numeric")
    r <- if(isNum) numeric(nz) else complex(nz)
    for(i in seq_len(nz)) {
	## zairy(zr, zi, id, kode, air, aii, nz, ierr)
	ri <- .Fortran("zbiry",
		       zr[i], zi[i],
		       id = deriv,
		       kode= as.integer(1L + as.logical(expon.scaled)),
					# 1 or 2, exactly as desired
		       bir = double(1),
		       bii = double(1),
		       nz   = integer(1),
		       ierr = integer(1), PACKAGE = "Bessel")
	if(ri$ierr) {
	    f.x <- sprintf("'zairy(%g %s %gi, deriv=%d)'", zr[i],
                           c("-","+")[1+(z >= 0)], abs(zi[i]), deriv)
	    if(ri$ierr == 3)
		warning(sprintf("%s large arguments -> precision loss\n	 (of at least half machine accuracy)", f.x))
	    else if(ri$ierr == 2) {
		if(getOption("verbose"))
                    message(sprintf("%s  -> overflow ; returning Inf\n", f.x))
                ri$bir <- ri$bii <- Inf
            }
	    else stop(sprintf("%s [Fortran] error ierr = %d", f.x, ri$ierr))
	}
	r[i] <- if(isNum) ri$bir else complex(real = ri$bir, imaginary = ri$bii)
    }
    r
} ## AiryB()
