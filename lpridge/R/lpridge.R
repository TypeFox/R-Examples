## lpridge	local polynomial regression fitting with ridging

lpridge <- function(x,y,bandwidth,
		 deriv = 0, n.out = 200,x.out = NULL, order = NULL,
		 ridge = NULL, weight = "epa", mnew = 100, var = FALSE)
{
    ## x		inputgrid
    ## y		data

    ## control and sort inputgrid and data
    n <- length(x)
    if (length(y) != n) stop("Input grid and data must have the same length.")
    sorvec <- sort.list(x)
    x <- x[sorvec]
    y <- y[sorvec]

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
	n.out <- length(x.out)
    }

    ## compute vector of bandwiths
    if (length(bandwidth) == 1)
	bandwidth <- rep(bandwidth,n.out)
    else if (length(bandwidth) != n.out)
	stop("Length of bandwith is not equal to length of output grid.")

    ## sort outputgrid and bandwidth
    sorvec <- sort.list(x.out)
    x.out <- x.out[sorvec]
    bandwidth <- bandwidth[sorvec]

    ## order	order of local polynomial approximation
    ## check order, deriv
    if (is.null(order)) order <- deriv+1
    if (order < 0) stop("Polynomial order is negative.")
    if (deriv < 0) stop("Order of derivative is negative.")
    if (deriv > order)
	stop("Order of derivative is larger than polynomial order.")

    ## ridge	ridge parameter
    if (is.null(ridge)) {
	ridge <- 5*sqrt(length(x)*mean(bandwidth)/diff(range(x)))*
	    mean(bandwidth)^(2*deriv)/((2*deriv+3)*(2*deriv+5))
	if (order == deriv) ridge <- 0
    }
    if (order == deriv & ridge > 0)
	stop("ridging is impossible for order==deriv.")

    ## weight	kernel weights
    if (weight == "epa") {
	kord <- 2
	wk <- c(1,0,-1)
    }
    else if (weight == "bi") {
	kord <- 4
	wk <- c(1,0,-2,0,1)
    }
    else if (weight == "tri") {
	kord <- 6
	wk <- c(1,0,-3,0,3,0,-1)
    }
    else if (is.numeric(weight)) {
	kord <- length(weight)-1
	wk <- weight
    }
    else stop("Error in weight.")

    ## mnew		force of restart

    ## var		switch for variance estimation
    var <- as.logical(var)

    ## internal parameters and arrays (see fortran routine)
    leng <- 10
    nmoms <- as.integer(length(x)/leng+1)
    imoms <- integer(nmoms)
    moms <- double(nmoms*4*(order+max(2,kord)+as.integer(var)))

    ## check internal limitations from fortran routine
    if (order > 10)
	stop("polynomial order exceeds 10.")
    if ((kord+order) > 12)
	stop("Order of kernel weights + polynomial order exceeds 12.")

    res <- .Fortran(lpridge_s,
		    x = as.double(x),
		    y = as.double(y),
		    as.integer(n),
		    bandwidth = as.double(bandwidth),
		    deriv = as.integer(deriv),
		    order = as.integer(order),
		    kord = as.integer(kord),
		    wk = as.double(wk),
		    x.out = as.double(x.out),
		    as.integer(n.out),
		    mnew = as.integer(mnew),
		    as.integer(imoms),
		    as.double(moms),
		    est = double(n.out),
		    as.integer(leng),
		    as.integer(nmoms),
		    var = as.integer(var),
		    est.var = double(n.out),
		    ridge = as.double(ridge),
		    nsins = integer(1))[c("bandwidth", "est","est.var", "nsins")]

    if (res$nsins > 0)
	warning(res$nsins,
		" singularity exceptions. Corresponding estimators set zero.")

    list(x = x, y = y, bandwidth = res$bandwidth, deriv = deriv, x.out = x.out,
	 order = order, ridge = ridge, weight = wk, mnew = mnew,
	 var = var, est = res$est, est.var = res$est.var)
}

