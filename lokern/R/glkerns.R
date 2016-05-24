### glkerns   kernel regression smoothing with bandwidth selection

## auxiliary, factored out of glkerns()
.glkerns <- function(x,y,x.out,nobs,n.out,deriv,korder,
		     hetero,is.rand,inputb,
		     m1,xl,xu,s,sig,bandwidth, trace.lev)
{
    ## calling fortran routine
    r <- .Fortran(glkern_s,			 # Fortran arg.names :
                  x = as.double(x),              # t
                  y = as.double(y),              # x
                  x.out = as.double(x.out),      # tt
		  est	= double(n.out),	 # y
		  nobs = as.integer(nobs),	 # n
                  n.out= as.integer(n.out),      # m
                  deriv = as.integer(deriv),     # nue
                  korder = as.integer(korder),   # kord
                  hetero = as.logical(hetero),   # hetero
                  is.rand= as.logical(is.rand),  # isrand
		  inputb = as.logical(inputb),	 # inputb
		  iter = as.integer(m1),# number of plug-in iterations on output
                  xl = as.double(xl),
                  xu = as.double(xu),
                  s = as.double(s),
                  sig = as.double(sig),
                  work1 = double((nobs+1)*5),	# wn [0:n, 5]
                  work2 = double(3 * m1),	# w1 [ m1, 3]
		  bandwidth = as.double(bandwidth)# = 19
		  , as.integer(trace.lev)
		  )[-c(1:2, 17:18, 20L)]	# all but (x,y), work*,..
    if(r$korder != korder)
	warning(gettextf("'korder' reset from %d to %d, internally",
			 korder, r$korder))
    if(r$iter < 0) r$iter <- NA_integer_
    r
}

glkerns <- function(x, y=NULL, deriv = 0,
                    n.out = 300, x.out = NULL, x.inOut = TRUE,
		    korder = deriv + 2, hetero = FALSE, is.rand = TRUE,
		    inputb = is.numeric(bandwidth) && bandwidth > 0,
		    m1 = 400, xl = NULL, xu = NULL, s = NULL, sig = NULL,
		    bandwidth = NULL, trace.lev = 0)
{
    ## control and sort input (x,y) - new: allowing only y
    xy <- xy.coords(x,y)
    x <- xy$x
    n <- length(x)
    if (n < 3) stop("must have n >= 3 observations")
    x.isInd <- !is.null(xy$xlab) && xy$xlab == "Index"
    isOrd <- x.isInd || !is.unsorted(x)
    if(isOrd)
        y <- xy$y
    else {
        ord <- sort.list(x)
        x <- x[ord]
	y <- y[ord]
    }

    ## compute/sort outputgrid 'x.out' (n.out : length of outputgrid)

    if (is.null(x.out)) {
	n.out <- as.integer(n.out)
	if(identical(x.inOut, FALSE)) {
	    x.out <- seq(x[1], x[n], length = n.out)
	}
	else { ## construct x.out containing x[] and more
	    if(identical(x.inOut, TRUE))
		seqXmethod <- "aim"	# cheaper than "interpolate"
	    else {			## x.inOut is a character
		if(!is.character(x.inOut))
		    stop("'x.inOut' must be logical or character")
		seqXmethod <- match.arg(x.inOut, eval(formals(seqXtend)$method))
		x.inOut <- TRUE
	    }
	    n.out <- length(x.out <- seqXtend(x, n.out, method = seqXmethod))
	    ind.x <- match(x, x.out)
	}
    }
    else {
	n.out <- length(x.out <- sort(x.out))
	ind.x <- match(x, x.out)## x[] matching x.out[]:
	## FIXME: approximate matching would be better: findInterval() etc
	if(!missing(x.inOut))
	    warning("'x.inOut' is disregarded as 'x.out' has been specified")
	x.inOut <- all(!is.na(ind.x))
	if(x.inOut)
	    seqXmethod <- NA_character_ # don't need seqXtend()
    }
    if(n.out == 0) stop("Must have 'n.out' >= 1")

    ## hetero	homo- or heteroszedasticity of error variables
    ## is.rand	random or non-random t-grid
    ## inputb	input bandwidth or estimation of plug-in bandwidth

    ## m1 : discretization for integral functional estimation
    if ((m1 <- as.integer(m1)) < 3)# was "10", but fortran has 3
        stop("number of discretizations 'm1' is too small")

    ## xl, xu: lower/upper bound for integral approximation and
    ##		variance estimation
    if (is.null(xl) || is.null(xu)) {
        xl <- 1
        xu <- 0
    }

    ## s	mid-point grid :
    s <- double(if(is.null(s) || length(s) != n+1)  n+1 else s)

    ## sig          input variance
    if (is.null(sig)) sig <- 0. #-> Fortran takes 0 = "compute default"

    inputb <- as.logical(inputb)
    if (is.null(bandwidth) || bandwidth < 0)
        bandwidth <- 0.
    else {
        bandwidth <- as.double(bandwidth[1])
        if (bandwidth == 0 && inputb)
            stop("bandwidth = 0 must have inputb = FALSE")
    }

    ## deriv          derivative of regression function to be estimated
    ## korder         kernel order
    if (deriv < 0 || deriv > 4)
        stop("Order of derivative 'deriv' must be in {0,1,..,4}.")
    if (deriv > 2 && !inputb)
        stop("Order of derivative must be <= 2  if (! inputb).")
    if (is.null(korder))
        korder <- deriv+2
    else if (korder > 6) {
        warning("Kernel order 'korder' must be <= 6; set to deriv + 2")
        korder <- deriv+2
    } else if (korder > 4 && !inputb) {
        warning("Kernel order must be <= 4 if(!inputb); set to deriv+2")
        korder <- deriv+2
    }

    xinL <- if(x.inOut) list(ind.x = ind.x, seqXmethod = seqXmethod)
    structure(c(xy[c("x","y")], # (x,y) possibly unsorted..
		.glkerns(x=x,y=y,x.out=x.out,nobs=n,n.out=n.out,deriv=deriv,
			 korder=korder,hetero=hetero,is.rand=is.rand,
			 inputb=inputb,m1=m1,xl=xl,xu=xu,
			 s=s,sig=sig,bandwidth=bandwidth,trace.lev=trace.lev),
		xinL,
		list(m1 = m1, isOrd = isOrd, ord = if(!isOrd) ord,
		     x.inOut = x.inOut, call = match.call())),
	      class = c("glkerns", "KernS"))
}
