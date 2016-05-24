qqevm <- function(object, nsim=1000, alpha=.050){
    # Want parameters as a matrix with one for for passing
    # through to family$rng etc.
    a <- t(object$coefficients)
    u <- object$threshold
    dat <- object$data$y

    qfun <- object$family$quant
    rfun <- object$family$rng

    ModPoints <- qfun(ppoints(dat), a, object)

    # If doing the envelope, simulate, sort and get the quantiles
    if (nsim > 0){
        n <- length(dat)
        sim <- matrix(rfun(nsim * n, a, object), ncol = nsim)
        sim <- apply(sim, 2, sort)
        # Get the simulated MSEs
        sim.mse <- apply(sim, 2, function(x, m) mean((x - m)^2), m = ModPoints)
        sim <- apply(sim, 1, quantile, prob = c(alpha/2, 1 - alpha/2))
    }
    else { sim <- NULL }

	# Get the quantile of the observed data MSE relative to the sims
    p <- mean(mean((sort(dat) - ModPoints )^2)  > sim.mse)

    res <- list(ModPoints=ModPoints, sim=sim, qfun=qfun, rfun=rfun, model=object,
                dat=dat, p=p)
    oldClass(res) <- 'qqevm'
    res
}

plot.qqevm <- function(x, xlab, ylab, main , plot = TRUE,
                       ylim = "auto", 
                       pch= 1, col =2 , cex=.75, linecol = 4 ,
                       intcol = 0, polycol = 15, ...){

    oldpar <- par(pty = "s"); on.exit(oldpar)	

    if (missing(xlab) || is.null(xlab)) { xlab <- "Model" }
    if (missing(ylab) || is.null(ylab)) { ylab <- "Empirical" }
    if (missing(main) || is.null(main)) { main <- "Quantile Plot" }

    limits <- if (is.null(x$sim)) NULL else range(x$sim, x$dat)

    plot(x$ModPoints, sort(x$dat),
         xlim = limits, ylim=limits,
         ylab = ylab, xlab = xlab, main = main,
         type = "n")
		# If doing the envelope, plot it before putting the data on 
    if (!is.null(x$sim)){
        if (polycol != 0){
            polygon(c(x$ModPoints,rev(x$ModPoints)),
                    c(x$sim[1, ], rev(x$sim[2, ])), col=polycol, border=NA)
        }
        lines(x$ModPoints, x$sim[1, ], col = intcol) 
        lines(x$ModPoints, x$sim[2, ], col = intcol)
    }
			
    # Add the diagonal reference line and the data
    abline(0, 1, col = linecol)
    points(x$ModPoints, sort(x$dat), pch = pch, col = col, cex=cex)
    box()

    invisible()
}

