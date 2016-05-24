trans <- function(object, ...) {
    UseMethod("trans", object)
} # end of 'trans' function.

datagrabber.extremalindex <- function(x, ...) {
    a <- attributes(x)
    # if(length(a$data.name) == 1) res <- c(get(a$data.name, ...))
    if(length(a$data.name) == 1) res <- c(a$data)
    else {
	look <- try(get(a$data.name[1], ...))
	if(class(look) == "try-error") look <- a$data
	nm <- colnames(look)
	if(!is.null(nm) & is.element(a$data.name[2], nm)) res <- look[,a$data.name[2]]
	else res <- look[,as.numeric(a$data.name[2])]
    }
    res <- a$na.action(res)
    return(res)
} # end of 'datagrabber.extremalindex' function.

datagrabber.declustered <- function(x, ...) {
    a <- attributes(x)

    if(length(a$data.name) == 1) {
	res <- c(try(get(a$data.name, ...), silent=TRUE))
	if(class(res) == "try-error") {
	    # res <- NULL
	    res <- a$data
	    return(res)
	}
    } else {
	look <- try(get(a$data.name[1], ...), silent=TRUE)
	if(class(look) == "try-error") {
	    # res <- NULL
	    res <- a$data
	    return(res)
	}
	nm <- colnames(look)
	if(!is.null(nm)) {
	    if(!is.element(a$data.name[2], nm)) dc <- as.numeric(a$data.name[2])
	    else dc <- a$data.name[2]

	    if(!is.element(a$data.name[3], nm)) dc <- c(dc, as.numeric(a$data.name[3]))
            else dc <- c(dc, a$data.name[3])
	    
	} else dc <- as.numeric(a$data.name[2:3])
	res <- cbind(c(look[,dc[1]]), c(look[,dc[2]]))
    }
    return(res)
} # end of 'datagrabber.declustered' function.

datagrabber.fevd <- function(x, response=TRUE, cov.data=TRUE, ...) {
    # Get response and any covariate data sets.
    in.data <- x$in.data
    if(response) {
	# if(!in.data) {
         #    if(!is.null(x$data.pointer)) y <- get(x$data.pointer, ...)
          #   else y <- x$x
	# } else y <- model.response(model.frame(x$x.fun, data=get(x$data.name[2], ...)))
	y <- x$x
    } else y <- NULL

    if(cov.data) {
        # if(!is.null(x$cov.data)) cdata <- x$cov.data
        # else if(x$data.name[2] == "") cdata <- NULL
        # else cdata <- get(x$cov.pointer, ...)
	cdata <- x$cov.data
    } else cdata <- NULL

    if(response & cov.data) out <- cbind(y, cdata)
    else if(response) out <- y
    else out <- cdata
    if(!is.null(out)) out <- do.call(x$na.action, list(out))
    return(out)
} # end of 'datagrabber.fevd' function.

trans.fevd <- function(object, ..., burn.in=499, return.all = FALSE) {

    x <- object
    n <- x$n

    meth <- tolower(x$method)

    y <- c(datagrabber(x)[,1])

    type <- tolower(x$type)
    des <- setup.design(x)

    if(is.element(meth, c("mle","gmle"))) p <- x$results$par
    else if(meth == "bayesian") {

	p2 <- x$results
	np <- dim(p2)[2] - 1
	if(np == 1) p2 <- matrix(p2[,1], ncol=1)
	else p2 <- p2[,1:np]
	if(burn.in != 0) p <- colMeans(p2[-(1:burn.in),])
	else p <- colMeans(p2)

    }

    pnames <- names(p)

    # Calculate the location parameters for each variable.
    if(is.element(type, c("gev", "pp", "gumbel", "weibull", "frechet"))) {
	X.loc <- des$X.loc
	nloc <- ncol(X.loc)
	loc <- rowSums(matrix(p[1:nloc], n, nloc, byrow=TRUE)*X.loc)
    } else nloc <- 0

    # Scale
    X.sc <- des$X.sc
    nsc <- ncol(X.sc)
    scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE)*X.sc)
    if(x$par.models$log.scale) scale <- exp(scale)

    # Shape
    if(!is.element(type, c("gumbel","exponential"))) {
	X.sh <- des$X.sh
 	nsh <- ncol(X.sh)
	shape <- rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], n, nsh, byrow=TRUE)*X.sh)
    } else shape <- 0

    if(is.element(type, c("pp", "gp", "exponential", "beta", "pareto"))) {
	u <- x$threshold
	eid <- y > u
	if(!return.all) {
	    scale <- scale[eid]
	    y <- y[eid]
	    if(length(u) > 1) u <- u[eid]
	    shape <- shape[eid]
	    if(type=="pp") scale <- scale + shape*(u - loc[eid])
	} else if(type=="pp") scale <- scale + shape*(u - loc)
    }

    if(is.element(type, c("gev", "weibull", "frechet"))) z <- -log(as.vector((1 + (shape * (y - loc))/scale)^(-1/shape)))
    else if(type=="gumbel") z <- as.vector((y - loc)/scale)
    else if(is.element(type, c("gp","beta","pareto"))) z <- -log(as.vector((1 + (shape * (y - u))/scale)^(-1/shape)))
    else if(type=="pp") z <- as.vector((1 + (shape * (y - u)/scale))^(-1/shape))
    else if(type=="exponential") z <- as.vector((y - u)/scale)
    return(z)

} # end of 'trans.fevd' function.

revtrans.evd <- function(z, threshold=NULL, location=NULL, scale, shape=NULL, type=c("GEV","GP","PP","Gumbel","Weibull","Frechet","Exponential","Beta","Pareto")) {
    type <- match.arg(type)
    if(is.element(type,c("Weibull","Frechet"))) type <- "GEV"
    else if(is.element(type,c("Beta","Pareto"))) type <- "GP"
    type <- tolower(type)
    if(is.element(type, c("gev","pp","gp"))) {
	zid <- shape==0
	if(any(zid)) shape[zid] <- 1e-10
    }
    if(type=="gumbel") out <- scale * z + location
    else if(type=="exponential") out <- scale * z + threshold
    else if(type=="gev") out <- location + (scale/shape) * (exp(shape * z) - 1)
    else if(type=="gp") out <- threshold + (scale/shape) * (exp(shape * z) - 1)
    else if(type=="pp") out <- threshold - (1 + shape * log(z)) * (scale + shape * (threshold - location))/shape
    return(out)
} # end of 'revtrans.evd' function.

blockmaxxer <- function(x, ...) {
    UseMethod("blockmaxxer", x)
} # end of 'blockmaxxer' function.

blockmaxxer.fevd <- function(x, ...) {

    if(x$type != "PP") stop("blockmaxxer: method for PP models only.")

    y <- datagrabber(x)

    blocks <- rep(1:x$span, each=x$npy)
    n2 <- length(blocks)
    if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
    else if(n2 > x$n) blocks <- blocks[1:x$n]

    res <- blockmaxxer(y, blocks = blocks)

    return(res)

} # end of 'blockmaxxer.fevd' function.

blockmaxxer.vector <- function(x, ..., blocks = NULL, blen = NULL, span = NULL) {

    if(is.null(blocks) && (is.null(blen) || is.null(span))) stop("blockmaxxer: must supply one of blocks or blen and span.")

    n <- length(x)
    if(is.null(blocks)) {
        blocks <- rep(1:span, each = blen)
        nb <- length(blocks)
        if(nb < n) blocks <- c(blocks, rep(blocks[nb], n - nb))
        else if(nb > n) blocks <- blocks[1:n]
    }

    res <- c(aggregate(x, by = list(blocks), max, ...)$x)
    return(res)
} # end of 'blockmaxxer.vector' function.

blockmaxxer.data.frame <- function(x, ..., which = 1, blocks = NULL, blen = NULL, span = NULL) {

    if(is.null(blocks) && (is.null(blen) || is.null(span))) stop("blockmaxxer: must supply one of blocks or blen and span.")

    xd <- dim(x)

    if(is.null(blocks)) {
	blocks <- rep(1:span, each = blen)
	nb <- length(blocks)
	if(nb < xd[1]) blocks <- c(blocks, rep(blocks[nb], xd[1] - nb))
	else if(nb > xd[1]) blocks <- blocks[1:xd[1]]
    }

    res <- c(aggregate(x[[which]], by = list(blocks), max, ...)$x)

    if(xd[2] > 1) {

        bl <- unique(blocks)
	n <- length(bl)

	o <- 1:length(x[[which]])
	ind <- numeric(0)

	tmp <- x[[which]]

	for(i in 1:n) {
	    id <- blocks == bl[i]
	    o2 <- o[id]
	    x2 <- tmp[id]
	    id2 <- x2 == max(x2, ...)
	    id2 <- id2

	    ind <- c(ind, min(o2[id2]))

	} # end of for 'i' loop.

	res <- x[ind,]
    }

    return(res)

} # end of 'blockmaxxer.data.frame' function.

blockmaxxer.matrix <- function(x, ..., which = 1, blocks = NULL, blen = NULL, span = NULL) {

    x <- as.data.frame(x)

    res <- blockmaxxer(x, ..., which = which, blocks = blocks, blen = blen, span = span)
    res <- as.matrix(res)

    return(res)

} # end of 'blockmaxxer.matrix' function.
