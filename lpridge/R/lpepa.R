## lpepa	local polynomials with Epanechnikov weights
##		for regression functions and derivatives

lpepa <- function(x, y, bandwidth,
                  deriv = 0, n.out = 200, x.out = NULL,
                  order = deriv+1, mnew = 100, var = FALSE)
{
    ## x		inputgrid
    ## y		data

    ## control and sort inputgrid and data
    n <- length(x)
    if (length(y) != n)
        stop("Input grid and data must have the same length.")
    sorvec <- sort.list(x)
    x <- as.double(x)[sorvec]
    y <- as.double(y)[sorvec]

    ## bandwidth	bandwidth for estimation

    ## deriv	derivative of regression function to be estimated

    ## n.out	length of outputgrid
    ## x.out	outputgrid

    ## compute and control outputgrid
    if (is.null(x.out)) {
        n.out <- as.integer(n.out)
        x.out <- seq(min(x),max(x),length = n.out)
    }
    else {
        x.out <- as.double(x.out)
        n.out <- length(x.out)
    }

    ## compute vector of bandwidths
    if (length(bandwidth) == 1)
	bandwidth <- rep(as.double(bandwidth), n.out)
    else {
	if (length(bandwidth) != n.out)
	    stop("Length of bandwidth is not equal to length of output grid.")
	storage.mode(bandwidth) <- "double"
    }

    ## sort outputgrid and bandwidth
    sorvec <- sort.list(x.out)
    x.out <- x.out[sorvec]
    bandwidth <- bandwidth[sorvec]

    ## order	order of local polynomial approximation
    ## check order, deriv
    if (order < 0) stop("Polynomial order is negative.")
    if (deriv > order)
        stop("Order of derivative is larger than polynomial order.")

    ## var		switch for variance estimation
    var <- as.logical(var)

    ## check internal limitations from fortran routine
    if (2 + order > 12)
        stop("Polynomial order exceeds 10.")

    ## internal parameters and arrays (see code in ../src/lpepa.f)
    leng <- 10
    nmoms <- as.integer(n/leng + 1)
    res <- .Fortran(lpepa_s,
                    x,
                    y,
                    as.integer(n),
		    bandwidth = bandwidth,
                    deriv = as.integer(deriv),
                    order = as.integer(order),
                    x.out = x.out,
                    as.integer(n.out),
                    as.integer(mnew),    ## mnew : force of restart
                    integer(nmoms), # imoms
                    double(nmoms*4*(2+order+as.integer(var))), # moms
                    est = double(n.out),
                    as.integer(leng),
                    nmoms,
                    var = as.integer(var),
                    est.var = double(n.out),
		    DUP = FALSE)[c("bandwidth", "est", "est.var")]

    list(x = x, y = y, bandwidth = res$bandwidth,
         deriv = deriv, x.out = x.out, order = order,
         mnew = mnew, var = var, est = res$est, est.var = res$est.var)
}

