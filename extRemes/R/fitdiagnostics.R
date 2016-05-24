lr.test <- function(x, y, alpha=0.05, df=1, ...) {

    if(class(x) == "fevd") {

        if(!is.element(x$method, c("MLE", "GMLE"))) stop("lr.test: fit method must be MLE or GMLE")
        l1 <- x$results$value
        df1 <- length(x$results$par)
        dname <- x$data.name[1]

    } else if(is.numeric(x) && length(x)==1) l1 <- x
    else if(is.numeric(x) && length(x)==2) l1 <- x[1]
    else stop("lr.test: invalid x argument.  Must be a single number, length two numeric vector, or an fevd object.")

    if(class(x) != "fevd") dname <- deparse(substitute(x))

    if(class(y) == "fevd") {

        if(!is.element(y$method, c("MLE", "GMLE"))) stop("lr.test: fit method must be MLE or GMLE")
        l2 <- y$results$value
        df2 <- length(y$results$par)
        dname <- c(dname, y$data.name[1])

    } else if(is.numeric(y) && length(y) == 1) l2 <- y
    else stop("lr.test: invalid y argument.  Must be a single number or an fevd object.")

    if((class(x) == "fevd") && class(y) == "fevd") {

        if(df2 < df1) return(lr.test(x=y, y=x, alpha=alpha, df=df, ...))
        else df <- df2 - df1

    }

    if(class(y) != "fevd") dname <- c(dname, deparse(substitute(y)))

    names(dname) <- c("x", "y")

    STATISTIC <- -2*(l2 - l1)
    names(STATISTIC) <- "Likelihood-ratio"
    CRITVAL <- qchisq(alpha, df=df, lower.tail=FALSE)
    PVAL <- pchisq(STATISTIC, df=df, lower.tail=FALSE)
    PARAMETER <- c(CRITVAL, alpha, df)
    names(PARAMETER) <- c("chi-square critical value", "alpha", "Degrees of Freedom")

    structure(list(statistic=STATISTIC, parameter=PARAMETER, alternative="greater", p.value=PVAL, method="Likelihood-ratio Test", data.name=dname), class="htest")

} # end of 'lr.test' function.

plot.fevd <- function(x, type = c("primary", "probprob", "qq", "qq2", "Zplot", "hist", "density", "rl", "trace"),
                    rperiods = c(2, 5, 10, 20, 50, 80, 100, 120, 200, 250, 300, 500, 800), a=0, hist.args=NULL, density.args=NULL, d = NULL, ...) {

    x2 <- x

    if(is.element(x$method, c("MLE", "GMLE"))) class(x2) <- "fevd.mle"
    else class(x2) <- paste("fevd.", tolower(x$method), sep="")

    UseMethod("plot", x2)

} # end of 'plot.fevd' function.

plot.fevd.lmoments <- function(x, type=c("primary", "probprob", "qq", "qq2", "Zplot", "hist", "density", "rl", "trace"),
                    rperiods=c(2, 5, 10, 20, 50, 80, 100, 120, 200, 250, 300, 500, 800), a=0, hist.args=NULL, density.args=NULL, d = NULL, ...) {

    # type <- tolower(type)
    type <- match.arg(type)
    type <- tolower(type)

    if(type=="trace") {
	cat("\n", "No trace plot available for L-moments estimation.\n")
	invisible()
    }

    args <- list(...)

    model <- x$type

    p <- x$results
    pnames <- names(p)

    if(is.element("location",pnames)) loc <- p["location"]
    else loc <- 0

    scale <- p["scale"]

    if(is.element("shape",pnames)) shape <- p["shape"]
    else shape <- 0

    if(!is.null(x$threshold)) u <- x$threshold
    else u <- 0

    y <- datagrabber(x, cov.data=FALSE)

    if(is.element(model, c("PP","GP","Exponential","Beta","Pareto"))) eid <- y > u
    else eid <- !logical(x$n)

    m <- sum(eid)

    if(is.element(type, c("primary","probprob","qq","rl"))) xp <- ppoints(m, a = a)

    if(type=="probprob") probprob.plot.evd(xp=xp, y=y, model=model, loc=loc, scale=scale, shape=shape, u=u, tform=FALSE, eid=eid, obj=x, npy=x$npy, ...)

    if(type=="primary") {

	if(model != "PP") par(mfrow=c(2,2), oma=c(0,0,2,0))
	else par(mfrow = c(3, 2), oma = c(0,0,2,0))

    }
    
    if(is.element(type,c("primary","qq"))) quantquant.plot.evd(x=x, xp=xp, y=y, u=u, loc=loc, scale=scale, shape=shape, tform=FALSE, eid=eid, model=model, type=type, ...)

    if(is.element(type,c("primary","qq2"))) quantquant2.plot.evd(x=x, y=y, eid=eid, model=model, type=type)

    if(type=="hist") histplot.evd(x=x, y=y, u=u, loc=loc, scale=scale, shape=shape, eid=eid, model=model, hist.args=hist.args, ...)

    if(is.element(type, c("primary","density"))) densplot.evd(x=x, y=y, u=u, loc=loc, scale=scale, shape=shape, eid=eid,
	model=model, tform=FALSE, density.args=density.args, type=type, ...)

    if(is.element(type, c("primary", "Zplot")) && model == "PP") {

	if(type == "primary") {

	    eeplot(x = x, type = "Zplot", main = "Zplot", d = d, ...)

	} else eeplot(x = x, type = "Zplot", d = d, ...)

    }    

    if(is.element(type, c("primary","rl"))) rlplot.evd(x=x, xp=xp, y=y, u=u, eid=eid, rperiods=rperiods, tform=FALSE, model=model, type=type, a=a, ...)

    if(type=="primary") mtext(deparse(x$call), line=0.5, outer=TRUE)

    invisible()
} # end of 'plot.fevd.lmoments' function.

plot.fevd.bayesian <- function(x, type=c("primary", "probprob", "qq", "qq2", "Zplot", "hist", "density", "rl", "trace"),
                    rperiods=c(2, 5, 10, 20, 50, 80, 100, 120, 200, 250, 300, 500, 800), a=0, hist.args=NULL, density.args=NULL,
		    burn.in=499, d = NULL, ...) {

    type <- match.arg(type)
    if(type != "Zplot") type <- tolower(type)

    model <- x$type

    p <- p2 <- x$results
    np <- ncol(p) - 1
    iters <- nrow(p)
    p <- p2 <- p[,1:np]
    pnames <- colnames(p)

    if(is.element("log.scale",pnames)) {

	id <- pnames == "log.scale"
	p[,id] <- exp(p[,id])
	pnames[id] <- "scale"
	colnames(p) <- pnames

    }

    p <- apply(p[-(1:burn.in), ], 2, mean, na.rm=TRUE)

    tform <- !is.fixedfevd(x)

    if(tform && is.element(model, c("PP", "GP", "Beta", "Exponential", "Pareto")) && is.element(type, c("hist", "density")))
	stop("plot.fevd: invalid type argument for this model.")

    if(!tform) {

	if(is.element(model, c("PP","GEV","Gumbel","Weibull","Frechet"))) {

	    loc <- p["location"]
	    nloc <- 1

	} else loc <- nloc <- 0

	scale <- p["scale"]
	nsc <- 1

	if(!is.element(model, c("Gumbel","Exponential"))) {

	    shape <- p["shape"]
	    nsh <- 1

	} else nsh <- shape <- 0

	ytrans <- NULL

    } else {

	ytrans <- trans(x)
	if(model == "PP") ytransPP <- trans(x, return.all = TRUE)

        designs <- setup.design(x)

        if(is.element(model, c("PP","GEV","Gumbel","Weibull","Frechet"))) {

	    X.loc <- designs$X.loc
	    nloc <- ncol(X.loc)
	    loc <- rowSums(matrix(p[1:nloc], x$n, nloc, byrow=TRUE) * X.loc)

        } else {

	    loc <- 0
	    nloc <- 0

        }

        X.sc <- designs$X.sc
        nsc <- ncol(X.sc)
        scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], x$n, nsc, byrow=TRUE) * X.sc)
        if(x$par.models$log.scale) scale <- exp(scale)

        if(!is.element(model, c("Gumbel","Exponential"))) {

	    X.sh <- designs$X.sh
	    nsh <- ncol(X.sh)
	    shape <- rowSums(matrix(p[(nloc+nsc+1):np], x$n, nsh, byrow=TRUE) * X.sh)

        } else {

    	    shape <- 0
    	    nsh <- 0

        }

	if(is.element(type, c("primary","return.level"))) if(missing(rperiods)) rperiods <- c(2, 20, 100)
	
    } # end of if '!tform' stmts.

    if(!is.null(x$threshold)) u <- x$threshold
    else u <- 0

    y <- datagrabber(x, cov.data=FALSE)


    if(model == "PP" && is.element(type, c("primary", "hist", "density"))) {


        blocks <- rep(1:round(x$span), each=round(x$npy))

        n2 <- length(blocks)
        if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
        else if(n2 > x$n) blocks <- blocks[1:x$n]

        ybm <- c(aggregate(y, by=list(blocks), max)$x)

    }

    if(is.element(model, c("PP","GP","Exponential","Beta","Pareto"))) eid <- y > u
    else eid <- !logical(x$n)

    m <- sum(eid)

    if(is.element(type, c("primary","probprob","qq","rl"))) xp <- ppoints(m, a = a)

    if(type=="probprob") probprob.plot.evd(xp=xp, y=y, model=model, loc=loc, scale=scale, shape=shape, u=u, tform=tform, eid=eid, obj=x, ytrans=ytrans, npy=x$npy, ...)

    if(type=="primary") {

	if(model != "PP") par(mfrow=c(2,2), oma=c(0,0,2,0))
	else if(!tform) par(mfrow = c(3,2), oma = c(0,0,2,0))
	else par(mfrow = c(2,2), oma = c(0,0,2,0))

    } else if(type=="trace") par(mfcol=c(2,np), oma=c(0,0,2,0))

    if(is.element(type,c("primary","qq"))) quantquant.plot.evd(x=x, xp=xp, y=y, u=u, loc=loc, scale=scale, shape=shape,
								    tform=tform, ytrans=ytrans, eid=eid, model=model, type=type, ...)

    if(is.element(type,c("primary","qq2"))) quantquant2.plot.evd(x=x, y=y, eid=eid, model=model, type=type)

    if(type=="hist") {

	if(model != "PP") histplot.evd(x=x, y=y, u=u, loc=loc, scale=scale, shape=shape, ytrans=ytrans,
		tform=tform, eid=eid, model=model, hist.args=hist.args, ...)
	else if(!tform) histplot.evd(x=x, y=ybm, u=u, loc=loc, scale=scale, shape=shape, ytrans=ytrans,
                tform=tform, eid=eid, model=model, hist.args=hist.args, ...)

    }

    if(is.element(type, c("primary", "Zplot")) && model == "PP") {

        if(type == "primary") {

            eeplot(x = x, type = "Zplot", main = "Z plot", d = d, ...)

        } else eeplot(x = x, type = type, d = d, ...)

    }

    if(type == "primary" && tform && is.element(model, c("PP", "GP", "Beta", "Exponential", "Pareto"))) { 

	## Eric -- At some point, should change ths to be a little cleaner, but ok for now.
	## Just do not want to plot the density plot if it is a non-stationary POT model.

    } else if(is.element(type, c("primary","density"))) {

        if(model != "PP") densplot.evd(x=x, y=y, u=u, loc=loc, scale=scale, shape=shape, tform=tform, eid=eid, ytrans=ytrans,
								    model=model, density.args=density.args, type=type, ...)
	else if(!tform) densplot.evd(x=x, y=y, u=u, loc=loc, scale=scale, shape=shape, tform=tform, eid=eid, ytrans=ybm,
                                                                    model=model, density.args=density.args, type=type, ...)
    }
    

    if(is.element(type, c("primary","rl"))) {

	# Eric 8/27/13 -- make sure for POT models that if the threshold varies, so do parameters, or else don't do rl plot.

	do.rlplot <- TRUE

	if(is.element(model, c("PP", "GP", "Exponential", "Beta", "Pareto"))) {
	
	    const.thresh <- check.constant(x$par.models$threshold)
	    const.loc <- check.constant(x$par.models$location)
	    const.scale <- check.constant(x$par.models$scale)
	    const.shape <- check.constant(x$par.models$shape)

	    if(!const.thresh && all(c(const.loc, const.scale, const.shape))) do.rlplot <- FALSE

	    if(!do.rlplot && type == "rl")
		stop("plot.fevd: invalid type argument for POT models with varying thresholds but constant parameters (are you sure you about this model choice?).")

	}

	if(do.rlplot) rlplot.evd(x=x, xp=xp, y=y, u=u, eid=eid, rperiods=rperiods, tform=tform, model=model, type=type, a=a, ...)

    }

    if(type=="trace") {

	pnames2 <- colnames(p2)
	id <- pnames2 == "log.scale"
	if(any(id)) p2[,id] <- exp(p2[,id])

	msg <- paste("Posterior Density\n", pnames, sep="")
	for(i in 1:np) {

	    if(is.null(density.args)) {

		if(i==1) plot(density(p2[-(1:burn.in),i]), main=msg[i], ...)
		else plot(density(p2[-(1:burn.in),i]), ylab="", main=msg[i], ...)

	    } else {

		yd <- do.call("density", c(list(x=p2[-(1:burn.in),i]), density.args))
		if(i==1) plot(yd, main=msg[i], ...)
		else plot(yd, ylab="", main=msg[i], ...)

	    }

	    xt <- 1:iters
	    if(i==1) plot(xt, p2[,i], type="l", xlab=pnames[i], ylab="MCMC trace", ...)
	    else plot(xt, p2[,i], type="l", xlab=pnames[i], ylab="", ...)
	    abline(h=p[i], v=burn.in, lty=2, col="gray")

	    # xt2 <- xt[(burn.in+1):iters]
	    # mt <- (cumsum(p2[(burn.in+1):iters,i])/xt2) - p[i]
	    # par(usr=c(range(xt), 1.5 * range(mt,finite=TRUE)))
	    # lines(xt2, mt, col="darkorange", lty=3)
	    # axis(4, at=pretty(mt), labels=pretty(mt), col="darkorange", col.ticks="darkorange")
	    # legend("topleft", legend=c("trace","posterior mean","cumsum diff"), col=c("black","gray","darkorange"), lty=1:3, bty="n")

	} # end of for 'i' loop.
    } # end of if 'trace' stmts.

    if(is.element(type, c("primary","trace"))) mtext(deparse(x$call), line=0.5, outer=TRUE)

    invisible()

} # end of 'plot.fevd.bayesian' function.

plot.fevd.mle <- function(x, type=c("primary", "probprob", "qq", "qq2", "Zplot", "hist", "density", "rl", "trace"),
		    rperiods=c(2, 5, 10, 20, 50, 80, 100, 120, 200, 250, 300, 500, 800), a=0, hist.args=NULL, density.args=NULL, period="year", 
		    prange=NULL, d = NULL, ...) {

    type <- match.arg(type)
    model <- x$type
    pars <- x$results$par
    args <- list(...)

    # Get response and any covariate data sets.
    ytmp <- datagrabber(x)
    if(x$data.name[2] != "") data <- ytmp[,-1]
    else data <- NULL
    y <- c(ytmp[,1])

    const.thresh <- check.constant(x$par.models$threshold)
    const.loc <- check.constant(x$par.models$location)
    const.scale <- check.constant(x$par.models$scale)
    const.shape <- check.constant(x$par.models$shape)

    # Eric -- 8/27/13 -- Don't allow rl plot if model is POT, threshold varies, but parameters are all constant.
    if(is.element(model, c("PP", "GP", "Exponential", "Beta", "Pareto")) && !const.thresh && all(c(const.loc, const.scale, const.shape)) && type == "rl")
	stop("plot.fevd: invalid type argument for POT models with varying thresholds but constant parameters (are you sure you about this model choice?).")

    tform <- !is.fixedfevd(x)
    if(tform) {
	ytrans <- trans(x)
	if(model == "PP" && is.element(type, c("primary", "hist", "density"))) ytransPP <- trans(x, return.all = TRUE)
    }

    if(!tform) {
	if(!is.element(model, c("GP","Beta","Pareto","Exponential"))) loc <- pars["location"]
        else loc <- NULL
        if(is.element("scale", names(pars))) scale <- pars["scale"]
        else scale <- exp(pars["log.scale"])
        if(!is.element(model, c("Gumbel","Exponential"))) shape <- pars["shape"]
        else shape <- 0
    } else if(is.element(type, c("primary", "rl"))) if(missing(rperiods)) rperiods <- c(2, 20, 100)
    # end of if else '!tform' stmts.

    npy <- x$npy

    if(is.element(model, c("PP","GP","Beta","Pareto","Exponential"))) {
	u <- x$threshold
	eid <- y > u
	lam <- mean(eid)
    } else u <- lam <- NULL

    if(is.element(type,c("primary","probprob","qq", "rl"))) {

	if(is.element(model, c("GP","PP","Beta","Pareto","Exponential"))) {

            if(!tform) n <- sum(y > x$threshold)
            else n <- sum(!is.na(ytrans) & !is.nan(ytrans))

        } else n <- x$n

	xp <- ppoints(n = n, a = a)
    }

    if(type=="primary") {

	# if(is.element(model, c("GP", "Exponential", "Beta", "Pareto")) && tform) par(mfrow = c(1, 2), oma = c(0,0,2,0))
	if(model != "PP") par(mfrow=c(2,2), oma=c(0,0,2,0))
	else if(!tform) par(mfrow = c(3, 2), oma = c(0,0,2,0))
	else par(mfrow = c(2,2), oma = c(0,0,2,0))

    }

    if(type=="probprob") {
	if(is.null(args$main))  m1 <- deparse(x$call)

        if(!tform) {

	    if(is.element(model, c("PP","GP","Beta","Pareto","Exponential"))) {

		yp <- pevd(sort(y[eid]), loc=loc, scale=scale, shape=shape, threshold=u, npy=npy, type=model)

	    } else yp <- pevd(sort(y), loc=loc, scale=scale, shape=shape, threshold=u, npy=npy, type=model)

            if(is.null(args$main)) plot(xp, yp, main=m1, xlab="Empirical Probabilities", ylab="Model Probabilities", ...)
	    else  plot(xp, yp, xlab="Empirical Probabilities", ylab="Model Probabilities", ...)

	    abline(0,1)

	} else {
	    # yp <-  pevd(sort(ytrans), loc=0, scale=1, shape=0, threshold=0, npy=npy, type=model)
	    if(is.element(model, c("GEV","Gumbel","Weibull","Frechet"))) yp <- exp(-exp(-sort(ytrans)))
	    else if(model=="PP") yp <- sort(ytrans)
	    else yp <- 1 - exp(-sort(ytrans))
	    if(is.null(args$main)) {
		plot(xp, yp, main=m1, xlab="Residual Empirical Probabilities", ylab="Residual Model Probabilities", ...)
            } else  plot(xp, yp, xlab="Residual Empirical Probabilities", ylab="Residual Model Probabilities", ...)
	    abline(0,1)
	}
    } # end of 'probprob' stmts.

    if(is.element(type,c("primary","qq"))) {

	if(!tform) {

	    if(is.element(model, c("Weibull","Frechet"))) mod2 <- "GEV"
	    else if(is.element(model, c("Exponential","Beta","Pareto"))) mod2 <- "GP"
	    else mod2 <- model

	    if(!is.element(model,c("PP","GP","Beta","Pareto","Exponential"))) yq <- qevd(xp, loc=loc, scale=scale, shape=shape, type=mod2)
	    else yq <- qevd(xp, threshold=u, loc=loc, scale=scale, shape=shape, type=mod2)

	    if(is.null(args$main)) {

		if(type=="primary") m2 <- ""
		else m2 <- deparse(x$call)
	        if(is.element(model, c("GEV","Weibull","Frechet","Gumbel"))) plot(yq, sort(y), xlab="Model Quantiles", ylab="Empirical Quantiles", main=m2)
		else plot(yq, sort(y[eid]), xlab="Model Quantiles", ylab="Empirical Quantiles", main=m2, ...)

	    } else {

		if(is.element(model, c("GEV","Weibull","Frechet","Gumbel"))) plot(yq, sort(y), xlab="Model Quantiles", ylab="Empirical Quantiles", ...)
		else plot(yq, sort(y[eid]), xlab="Model Quantiles", ylab="Empirical Quantiles", ...)

	    }

	} else {

	    if(is.null(args$main)) {

		if(type=="primary") m2 <- ""
		else m2 <- deparse(x$call)

		if(is.element(model, c("GEV","Weibull","Frechet"))) m2 <- paste(m2, "(Gumbel Scale)", sep="\n")
		else if(is.element(model, c("PP", "GP", "Beta", "Pareto"))) m2 <-  paste(m2, "Exponential Scale", sep="\n")

	        if(is.element(model, c("GEV","Weibull","Gumbel","Frechet"))) plot(-log(-log(sort(xp))), sort(ytrans), main=m2,
										    xlab="(Standardized) Model Quantiles", ylab="Empirical Residual Quantiles", ...)
		else if(is.element(model, c("GP","Beta","Exponential","Pareto"))) plot(-log(1 - xp), sort(ytrans), main=m2,
										    xlab="(Standardized) Residual Quantiles", ylab="Empirical Residual Quantiles", ...)
		else if(model=="PP") plot(-log(1 - xp), sort(-log(ytrans)), main=m2, xlab="(Standardized) Residual Quantiles", ylab="Empirical Residual Quantiles", ...)
	    } else {
		if(is.element(model, c("GEV","Weibull","Gumbel","Frechet"))) plot(-log(-log(sort(xp))), sort(ytrans), xlab="Model", ylab="Empirical", ...)
		else if(is.element(model, c("GP","Beta","Exponential","Pareto"))) plot(-log(1 - xp), sort(ytrans), xlab="Model", ylab="Empirical", ...)
		else if(model=="PP") plot(-log(1 - xp), sort(-log(ytrans)), xlab="Model", ylab="Empirical", ...)
	    }
	} # end of if else 'tform' stmts.
	abline(0,1)
    } # end of if 'qq' stmts.

    if(is.element(type, c("primary","qq2"))) {

      if(is.fixedfevd(x) && !is.null(x$blocks)) {

        z <- rextRemes(x, round(x$blocks$nBlocks * x$npy)) # CJP: for stationary, it generates based on number of obs, so need to tell it how many obs there would be

      } else {

        z <- rextRemes(x)  # for non-stationary, it generates based on number of exceedances, so this is fine

      }
	if(!is.element(model, c("PP","GP","Beta","Exponential","Pareto"))) {

	    yQQ <- y
	    if(is.null(args$xlab)) xl <- paste(x$data.name[1], " Empirical Quantiles", sep="")
	    else xl <- args$xlab

	} else {

	    yQQ <- y[eid]
	    if(is.null(args$xlab)) {

		if(length(x$threshold)==1) xl <- paste(x$data.name[1], "( > ", x$threshold, ") Empirical Quantiles", sep="")
		else xl <- paste(x$data.name[1], "( > threshold) Empirical Quantiles", sep="")

	    } else xl <- args$xlab
	}

	if(!is.null(args$main)) {

	    if(type=="primary") mQQ <- ""
	    else mQQ <- deparse(x$call)
	    qqplot(yQQ, z, main=mQQ, xlab=xl, ylab="Quantiles from Model Simulated Data")

	} else {

	    qqplot(yQQ, z, xlab=xl, ylab="Quantiles from Model Simulated Data")

	}

    } # end of if do second qq-plot using simulated data from the model stmts.

    if(!(type == "primary" && model == "PP" && tform)) {

        if(is.element(type, c("primary","density")) && is.null(x$blocks)) {
    
     # CJP -- this needs the full data vector (or perhaps more extensive rewriting to work with blocks)
    
          if(!tform) {
    
            if(!is.element(model, c("PP","GP","Beta","Exponential","Pareto"))) yd <- y
            else if(model != "PP") yd <- y[eid] - x$threshold
            else {
              if(x$span %% 1 != 0 || x$npy %% 1 != 0)
                warning("plot.fevd.mle: span or npy not integers; determination of max in each block may be substantially in error when there are many blocks.")
              blocks <- rep(1:round(x$span), each=round(x$npy))
    
        # CJP2: hard to deal exactly with inhomogeneity in number of values in each block;
        # hopefully rounding here is reasonable in most cases but with many blocks,
        # the indexing could be substantially off.
    
              n2 <- length(blocks)
              if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
              else if(n2 > x$n) blocks <- blocks[1:x$n]
    
              yd <- c(aggregate(y, by=list(blocks), max)$x)
    
            }
   

            if(is.null(density.args)) yd <- density(yd)
            else yd <- do.call("density", c(list(yd), density.args))

	    if(is.element(type, c("primary","density"))) xd <- seq(min(yd$x, na.rm=TRUE), max(yd$x, na.rm=TRUE),,100)
            else xd <- seq(min(yh, na.rm=TRUE), max(yh, na.rm=TRUE),,100)
            if(is.element(model, c("PP", "Gumbel", "Weibull", "Frechet"))) mod2 <- "GEV"
            else if(is.element(model, c("Pareto", "Frechet", "Beta", "Exponential"))) mod2 <- "GP"
            else mod2 <- model

	    if(tform) {

                mu <- 0
                sig <- 1
                xi <- 0
                u2 <- 0

            } else {

                mu <- loc
                sig <- scale
                xi <- shape
                u2 <- u[1]
                # if(model=="PP") sig <- sig + xi * (u2 - loc)

            }

            yd2 <- devd(xd, loc=mu, scale=sig, shape=xi, threshold=u2, type=mod2)
            # lines(xd, yd2, lty=2, col="blue", lwd=1.5)
            
            if(is.null(args$ylim)) {

              yld <- range(c(yd$y, yd2), finite=TRUE)
              yld[1] <- min(yld[1], 0)
              # yld[2] <- max(yld[2]+0.5, 1)

            }
            
            if(is.null(args$main)) {

              if(type=="primary") m3 <- ""
              else m3 <- deparse(x$call)

              if(is.null(args$ylim)) plot(yd, main=m3, ylim=yld, ...)
              else plot(yd, main=m3, ...)

            } else {

              if(is.null(args$ylim)) plot(yd, ylim=yld, ...)
              else plot(yd, ...)

            }

          } else {
    
    	if(type == "density" && is.element(model, c("PP", "GP", "Exponential", "Pareto", "Beta")))
    	    stop("plot.fevd: invalid type argument for this model.")
    
            if(!is.element(model, c("PP", "GP", "Exponential", "Pareto", "Beta"))) {
    # 	if(model == "PP") { # Eric 8/14/13
    # 	    if(x$span %% 1 != 0 || x$npy %% 1 != 0)
    #             warning("plot.fevd.mle: span or npy not integers; determination of max in each block may be substantially in error when there are many blocks.")
    #           blocks <- rep(1:round(x$span), each=round(x$npy))
    # 
    #     # CJP2: hard to deal exactly with inhomogeneity in number of values in each block;
    #     # hopefully rounding here is reasonable in most cases but with many blocks,
    #     # the indexing could be substantially off
    # 
    # 
    #           n2 <- length(blocks)
    #           if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
    #           else if(n2 > x$n) blocks <- blocks[1:x$n]
    # 
    # 	    ytmp <- c(aggregate(y, by = list(blocks), max)$x)
    # 
    #           mbid <- logical(x$n)
    # 
    #           for(i in 1:(length(unique(blocks)))) {
    # 
    #               tmpind <- blocks == blocks[i]
    #               tmpind2 <- y == ytmp[i]
    #               tmpind2[is.na(tmpind2)] <- FALSE
    # 
    #               finind <- tmpind & tmpind2
    # 
    #               if(sum(finind) > 1) {
    # 
    #                   numsind <- (1:x$n)[finind]
    #                   numsind <- numsind[1]
    #                   finind[-numsind] <- FALSE
    # 
    #               }
    # 
    #               mbid[finind] <- TRUE
    # 
    #           } # end of for 'i' loop.
    # 
    #           ytrans2 <- ytransPP[mbid]
    #           # ytrans2 <- c(aggregate(ytransPP, by=list(blocks), max)$x)
    # 
    #         } else ytrans2 <- ytrans
   
	    
            if(is.null(density.args)) yd <- density(ytrans)
            else yd <- do.call("density", c(list(ytrans), density.args))
            
            if(is.null(args$ylim)) {
              yld <- range(yd$y, finite=TRUE)
              yld[1] <- min(yld[1], 0)
              # yld[2] <- max(yld[2]+0.5, 1)
            }
            
            if(is.null(args$main)) {
    
              if(type=="primary") m3 <- "Transformed Data"
              else m3 <- paste(deparse(x$call), "Transformed Data", sep="\n")
              if(is.null(args$ylim)) plot(yd, main=m3, ylim=yld, ...)
              else plot(yd, main=m3, ...)
    
            } else {
    
              if(is.null(args$ylim)) plot(yd, ylim=yld, ...)
              else plot(yd, ...)
    
            }
    
    	} # end of if model is non-stationary GEV stmts.
          } # end of if else '!tform' stmts.
    
        } # end of if 'primary plots or density' stmts.
    }

    if(type=="hist" && is.null(x$blocks)) {  # CJP2: I previously forgot to exclude this when there are blocks

      if(!tform) {

        if(is.null(args$main)) {

          if(model != "PP") m4 <- paste(deparse(x$call), "Histogram", sep="\n")
          else m4 <- paste(deparse(x$call), "\n", paste("Histogram (", x$period.basis, " maxima)", sep=""))

	    }

	    if(is.element(model, c("GP", "Beta", "Exponential", "Pareto"))) yh <- y[eid] - x$threshold
	    else if(model != "PP") yh <- y
	    else {
		blocks <- rep(1:x$span, each=x$npy)
                n2 <- length(blocks)
                if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
                else if(n2 > x$n) blocks <- blocks[1:x$n]
                yh <- c(aggregate(y, by=list(blocks), max)$x)
	    }

	} else {

	    if(is.element(model, c("PP", "GP", "Exponential", "Beta", "Pareto"))) stop("plot.fevd: invalid type argument for this model.")
	    if(is.null(args$main)) m4 <- paste(deparse(x$call), "Histogram of Transformed Data", sep="\n")
	    # if(model != "PP") yh <- ytrans
	    # else {

	# 	blocks <- rep(1:x$span, each=x$npy)
         #        n2 <- length(blocks)
          #       if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
           #      else if(n2 > x$n) blocks <- blocks[1:x$n]
            #     yh <- c(aggregate(ytransPP, by=list(blocks), max)$x)

	    # } # end of if 'model != PP' stmts.

	}

	if(is.null(hist.args)) {
	    if(!is.null(args$ylim)) {
 	        if(is.null(args$main)) {
	            if(is.null(args$col)) hist(yh, col="darkblue", freq=FALSE, breaks="FD", xlab=x$data.name[1], main=m4, ...)
	            else  hist(yh, freq=FALSE, breaks="FD", main=m4, xlab=x$data.name[1], ...)
	        } else {
	            if(is.null(args$col)) hist(yh, col="darkblue", freq=FALSE, breaks="FD", xlab=x$data.name[1], ...)
                    else  hist(yh, freq=FALSE, breaks="FD", xlab=x$data.name[1], ...)
	        }
	    } else {
		if(is.null(args$main)) {
                    if(is.null(args$col)) hist(yh, col="darkblue", freq=FALSE, breaks="FD", xlab=x$data.name[1], main=m4, ylim=c(0,1.5), ...)
                    else  hist(yh, freq=FALSE, breaks="FD", main=m4, xlab=x$data.name[1], ylim=c(0,1.5), ...)
                } else {
                    if(is.null(args$col)) hist(yh, col="darkblue", freq=FALSE, breaks="FD", xlab=x$data.name[1], ylim=c(0,1.5), ...)
                    else  hist(yh, freq=FALSE, breaks="FD", xlab=x$data.name[1], ylim=c(0,1.5), ...)
                }
	    }
	} else do.call("hist", c(list(yh), hist.args))
    } # end of if 'type = hist' stmts.

    if(!(tform && is.element(model, c("PP", "GP", "Beta", "Exponential", "Pareto")))) {

        if(is.element(type, c("primary", "density", "hist")) && is.null(x$blocks)) {
    
            # CJP: this needs the full data vector (or perhaps more extensive rewriting to work with blocks)
    
    	    if(is.element(type, c("primary","density"))) xd <- seq(min(yd$x, na.rm=TRUE), max(yd$x, na.rm=TRUE),,100)
    	    else xd <- seq(min(yh, na.rm=TRUE), max(yh, na.rm=TRUE),,100)
    	    if(is.element(model, c("PP", "Gumbel", "Weibull", "Frechet"))) mod2 <- "GEV"
    	    else if(is.element(model, c("Pareto", "Frechet", "Beta", "Exponential"))) mod2 <- "GP"
    	    else mod2 <- model
    
    	    if(tform) {
    
    	        mu <- 0
    	        sig <- 1
    	        xi <- 0
    	        u2 <- 0
    
    	    } else {
    
    	        mu <- loc
    	        sig <- scale
    	        xi <- shape
    	        u2 <- u[1]
    	        # if(model=="PP") sig <- sig + xi * (u2 - loc)
    
    	    }

    	    yd2 <- devd(xd, loc=mu, scale=sig, shape=xi, threshold=u2, type=mod2) 
    	    lines(xd, yd2, lty=2, col="blue", lwd=1.5)
    
    	    if(model != "PP") {
    
    	        if(type=="hist" || !is.null(density.args)) legend("topright", legend="Modeled Density", col="blue", lty=2, lwd=1.5, bty="n")
    	        else if(is.null(density.args)) legend("topright", legend=c("Empirical", "Modeled"), col=c("black","blue"), lty=c(1,2), lwd=c(1, 1.5), bty="n")
    
    	    } else {
    
    	        if(type=="hist" || !is.null(density.args)) legend("topright", legend="Modeled Density", col="blue", lty=2, lwd=1.5, bty="n")
                else if(is.null(density.args)) legend("topright", legend=c(paste("Empirical (", x$period.basis, " maxima)", sep=""), "Modeled"),
    						    col=c("black","blue"), lty=c(1,2), lwd=c(1, 1.5), bty="n")
    
    	    }

        } # end of if add density lines stmts.

    } # end of making sure that it is not a non-stationary POT model stmts.


    # Z and W plots for all PP models.
    if(is.element(type, c("primary", "Zplot")) && model == "PP") {

        if(type == "primary") {

            eeplot(x = x, type = "Zplot", main = "Z plot", d = d, ...)

        } else eeplot(x = x, type = type, d = d, ...)

    }

    if(type == "Zplot" && model != "PP") stop("plot.fevd.mle: invalid type argument for this model.")


    if(is.element(type, c("primary", "rl")) && is.null(x$blocks)) {

        # Eric 8/27/13 -- do not plot return levels if model is POT, threshold varies and parameters do not vary.
        if(!(is.element(model, c("PP", "GP", "Exponential", "Beta", "Pareto")) && !const.thresh && all(c(const.loc, const.scale, const.shape)))) {

    # CJP: this needs the full data vector
    # (or perhaps more extensive rewriting to work with blocks,
    # but would be hard to deal with blocks with no exceedances)

	if(is.null(args$main)) {

	    if(type=="primary") m5 <- ""
	    else m5 <- deparse(x$call)

	}

	if(model=="PP") {

            # if(!tform) mod2 <- "GEV"
	    # else mod2 <- "GP"
	    mod2 <- "GEV"

            if(is.null(args$main)) {
		# if(!tform) {
		    if(type=="rl") m5 <- paste(m5, "Return Levels based on approx. equivalent GEV", sep="\n")
		    else m5 <- paste("Return Levels based on approx.", "equivalent GEV", sep="\n")
		# } else {
		 #    if(type=="rl") m5 <- paste(m5, "Return Levels based on approx. equivalent GP df", sep="\n")
                  #   else m5 <- paste("Return Levels based on approx.", "equivalent GP df", sep="\n")
		# }
	    }

	    blocks <- rep(1:x$span, each=x$npy)
	    n2 <- length(blocks)
	    if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
	    else if(n2 > x$n) blocks <- blocks[1:x$n]
	    yEmp <- c(aggregate(y, by=list(blocks), max)$x)

        } else mod2 <- model

	if(is.null(x$units)) ylb <- "Return Level"
        else ylb <- paste("Return Level (", x$units, ")", sep="")

	if(!tform) {

	    yrl <- rlevd(rperiods, loc=loc, scale=scale, shape=shape, threshold=u, type=mod2, npy=npy, rate=lam)
	    bds <- ci(x, return.period=rperiods)
	    if(is.null(args$ylim)) yl <- range(c(bds), finite=TRUE)
	    
	    if(is.element(model, c("PP", "GEV", "Gumbel", "Weibull", "Frechet"))) xrl <- -1/(log(1 - 1/rperiods))
	    else xrl <- rperiods

	    # if(is.null(x$units)) ylb <- "Return Level"
            # else ylb <- paste("Return Level (", x$units, ")", sep="")

	    xlb <- paste("Return Period (", x$period.basis, "s)", sep="")

	    if(is.null(args$main)) {
		if(!is.null(args$ylim)) plot(xrl, yrl, type="l", log="x", xlab=xlb, ylab=ylb, main=m5, ...)
		else plot(xrl, yrl, type="l", log="x", ylim=yl, xlab=xlb, ylab=ylb, main=m5, ...)
	    } else {
		if(!is.null(args$ylim)) plot(xrl, yrl, type="l", log="x", xlab=xlb, ylab=ylb, ...)
		else plot(xrl, yrl, type="l", log="x", ylim=yl, xlab=xlb, ylab=ylb, ...)
	    }

	    lines(xrl, bds[,1], col="gray", lty=2, lwd=2)
	    lines(xrl, bds[,3], col="gray", lty=2, lwd=2)

	    if(is.element(model, c("GEV", "Gumbel", "Weibull", "Frechet"))) points(-1/log(xp), sort(y))
	    else if(is.element(model, c("GP", "Beta", "Pareto", "Exponential"))) {
		n2 <- x$n
		if(is.null(a)) xp2 <- ppoints(n2)
		else xp2 <- ppoints(n2, a=a)
		sdat <- sort(y)
		points(-1/log(xp2)[sdat > u]/npy, sdat[sdat > u])
	    } else if(model == "PP") {
		if(is.null(a)) xp2 <- ppoints(length(yEmp))
		else xp2 <- ppoints(length(yEmp), a=a)
		points(-1/log(xp2), sort(yEmp))
	    }
	} else {

	    np <- length(rperiods)
	    effrl <- matrix(NA, length(y), np)
	    for(i in 1:np) effrl[,i] <- erlevd(x = x, period = rperiods[i])

	    if(is.null(args$ylim)) {

		yl <- range(c(c(y), c(effrl)), finite=TRUE)
		yl[2] <- yl[2] + sign(yl[2]) * log2(abs(yl[2]))

		if(is.null(args$main)) {

                    # if(!is.element(model, c("GP", "Beta", "Exponential", "Pareto")))
		    plot(y, type="l", xlab="index", ylab=ylb, main=m5, ylim=yl, ...)
		    # else plot(y[eid], type="l", xlab="index", ylab=ylb, main=m5, ylim=yl, ...)

                } else {

                    # if(!is.element(model, c("GP", "Beta","Exponential","Pareto")))
		    plot(y, type="l", xlab="index", ylab=ylb, ...)
                    # else plot(y[eid], type="l", xlab="index", ylab=ylb, ylim=yl, ...)

                }

	    } else {

	        if(is.null(args$main)) {

		    # if(!is.element(model, c("GP", "Beta","Exponential","Pareto")))
		    plot(y, type="l", xlab="index", ylab=ylb, main=m5, ...)
		    # else plot(y[eid], type="l", xlab="index", ylab=ylb, main=m5, ...)

	        } else {

		    # if(!is.element(model, c("GP", "Beta", "Exponential", "Pareto")))
		    plot(y, type="l", xlab="index", ylab=ylb, ...)
	    	    # else plot(y[eid], type="l", xlab="index", ylab=ylb, ...)

	        }

	    } # end of if else 'ylim passed via '...' stmts.

	    for(i in 1:np) lines(effrl[,i], lty=i, col=i+1)

	    if(is.element(model, c("PP", "GP", "Beta", "Exponential", "Pareto"))) {

                if(length(x$threshold) == 1) abline(h=x$threshold, col="darkorange", lwd=2)
                else lines(x$threshold, col="darkorange", lwd=2)

            }

	    if(!is.element(model, c("GEV", "Gumbel", "Weibull", "Frechet"))) {
		legend("topleft", legend=c(paste(rperiods, "-", period, " level", sep=""), "threshold"), lty=c(1:np,1), col=c(2:(np+1),"darkorange"), bg="white")
	    } else legend("topleft", legend=paste(rperiods, "-", period, " level", sep=""), lty=1:np, col=2:(np + 1), bg="white")

	} # end of if else '!tform' stmts.

	} # end of if model is not POT and threshold varies with non-varying parameters stmts.
    } # end of if 'make return level plots' stmts.

    if(type=="trace") {
	ntheta <- length(pars)
	if(is.null(prange)) {
	    theta.hat <- pars
	    if(check.constant(x$par.models$scale)) {
		phiU <- FALSE
		if(x$par.models$log.scale) theta.hat["log.scale"] <- exp(theta.hat["log.scale"])
		names(theta.hat)[names(theta.hat)=="log.scale"] <- "scale"
            } else phiU <- x$par.models$log.scale
	    prange <- matrix(NA, 2, ntheta)
	    tmp <- summary(x, silent=TRUE)
	    tmp <- rbind(tmp$par, tmp$se.theta)

	    if(is.matrix(tmp)) for(j in 1:ntheta) prange[,j] <- c(tmp[1,j] - 2*tmp[2,j], tmp[1,j] + 2*tmp[2,j])
	    else {
		tmp <- c(tmp)
		for(j in 1:ntheta) prange[,j] <- c(tmp[j] - 2*log2(abs(tmp[j])), tmp[j] + 2*log2(abs(tmp[j])))
	    }
	    if(any(tmp.id <- names(tmp)=="scale")) {
		if(!phiU) {
		    prange[1,tmp.id] <- max(prange[1,tmp.id], 1e-8)
		    if(prange[1,tmp.id] > prange[2, tmp.id]) prange[2, tmp.id] <- 2*prange[1, tmp.id]
		}
	    }
	} # end of if 'prange' is null or not stmts.

	hold <- list()
	par(mfrow=c(2,ntheta), oma=c(0,0,2,0))
	for(i in 1:ntheta) {  
	    if(!is.null(data)) hold[[i]] <- grlevdTracer(x=y, p=theta.hat, which.vary=i, p1.range=c(prange[,i]), threshold=u, threshold.fun=x$par.models$threshold,
						location.fun=x$par.models$location, scale.fun=x$par.models$scale, shape.fun=x$par.models$shape, data=data, phi=phiU, blocks=x$blocks, # CJP2 - changed u/threshold args
						type=model, npy=x$npy, na.action=x$na.action, par1.name=names(theta.hat)[i], plot=FALSE)
	    else hold[[i]] <- grlevdTracer(x=y, p=theta.hat, which.vary=i, p1.range=c(prange[,i]), threshold=u, threshold.fun=x$par.models$threshold, phi=phiU, blocks=x$blocks, # CJP2 - changed u/threshold arg calls
					        type=model, npy=x$npy, na.action=x$na.action, par1.name=names(theta.hat)[i], plot=FALSE)
	    if(i==1) plot(hold[[i]], type="likelihood", main="", xlab="")
	    else plot(hold[[i]], type="likelihood", main="", xlab="", ylab="")
	    abline(v=theta.hat[i], lty=2)
	} # end of for 'i' loop.

	for(i in 1:ntheta) {
	    if(i==1) plot(hold[[i]], type="gradient", main="")
	    else plot(hold[[i]], type="gradient", main="", ylab="")
	} # end of second for 'i' loop.
    } # end of if 'type == trace' stmts.

    if(is.element(type, c("primary","trace"))) mtext(deparse(x$call), line=0.5, outer=TRUE)
    if(type=="trace") invisible(hold)
    else invisible()

} # end of 'plot.fevd.mle' function.

erlevd <- function(x, period=100) {

  # CJP: this function now returns results per block rather than per observation when blocks is part of 'x'
    type <- x$type
    y <- datagrabber(x)
    if(x$data.name[2] != "") {
        data <- y[,-1]
        y <- y[,1]
    } else data <- NULL

    pars <- findpars(x, use.blocks = TRUE)  # CJP2: when x$blocks not NULL, findpars will now return pars with one row per block not per obs.

    if(is.null(pars$location)) pars$location <- rep(0, x$n)

    if(is.null(x$threshold)) u <- rep(0, x$n)
    else if(length(x$threshold)==1) u <- rep(x$threshold, x$n)
    else u <- x$threshold

    if(!is.null(x$threshold)) lam <- mean(y > u)
      else lam <- rep(1, x$n) # Eric -- changed 'len' to 'x$n'

    if(!is.null(x$blocks)) {  # CJP2
      u <- x$blocks$threshold
      if(is.null(pars$location)) pars$location <- rep(0, x$nBlocks) # this should never happen since when have blocks is only for PP model...
    }

    # CJP2: commented this out
    # if(type=="PP") {
    #   scale <- pars$scale
    #   scale <- scale + pars$shape * (u - pars$location)
    #} else scale <- pars$scale

    theta <- cbind(pars$location, pars$scale, pars$shape, u)  # CJP2: now 'pars$scale' not 'scale'
    # if(is.element(type,c("PP","GP","Beta","Exponential","Pareto"))) theta <- theta[y > u,]

    if(type=="PP") type2 <- "GEV"  # CJP2
    else type2 <- type

    ifun <- function(theta, pd, npy, rate) rlevd(period=pd, loc=theta[1], scale=theta[2], shape=theta[3], threshold=theta[4], type=type2, npy=npy, rate=rate)

    out <- apply(theta, 1, ifun, pd=period, npy=x$npy, rate=lam)

    return(out)
} # end of 'erlevd' function.


findpars <- function(x, ...) {
    UseMethod("findpars", x)
} # end of 'findpars' function.

findpars.fevd <- function(x, ...) {
    meth <- tolower(x$method)
    if(meth=="gmle") newcl <- "fevd.mle"
    else newcl <- paste("fevd.", meth, sep="")
    class(x) <- newcl
    UseMethod("findpars",x)
} # end of 'findpars.fevd' function.

findpars.fevd.lmoments <- function(x, ...) {
    return(x$results)
} # end of 'findpars.fevd.lmoments' function.

findpars.fevd.bayesian <- function(x, burn.in=499, FUN="mean", use.blocks=FALSE, ..., qcov = NULL) {

    model <- x$type

    tform <- !is.fixedfevd(x)

    f <- match.fun(FUN)

    p <- p2 <- x$results
    np <- ncol(p) - 1
    if(burn.in > 0) p <- p[-(1:burn.in),1:np]

    pnames <- colnames(p)

    if(is.element("log.scale", pnames)) {
        id <- pnames == "log.scale"
        p[,id] <- exp(p[,id])
        pnames[id] <- "scale"
        colnames(p) <- pnames
    }

    if(FUN=="mean") p <- colMeans(p)
    else p <- apply(p, 2, f)

    if(!tform) {

        if(is.element(model, c("PP","GEV","Gumbel","Weibull","Frechet"))) {
            loc <- p["location"]
            nloc <- 1
        } else loc <- NULL

        scale <- p["scale"]

        if(!is.element(model, c("Gumbel","Exponential"))) shape <- p["shape"]
        else shape <- NULL
 
       ytrans <- NULL

    } else {

	# Eric -- 8/8/13 -- I don't think this was used.
        # ytrans <- trans(x)

	if(is.null(qcov)) { # Eric 8/8/13 -- Attempting to allow for 'qcov' argument.

            if(is.null(x$blocks) || !use.blocks) { # CJP2
              designs <- setup.design(x)
          
              if(is.element(model, c("PP","GEV","Gumbel","Weibull","Frechet"))) {

                X.loc <- designs$X.loc
	        nloc <- ncol(X.loc)
                loc <- rowSums(matrix(p[1:nloc], x$n, nloc, byrow=TRUE) * X.loc)

              } else loc <- NULL
          
              X.sc <- designs$X.sc
              nsc <- ncol(X.sc)
              scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], x$n, nsc, byrow=TRUE) * X.sc)
              if(x$par.models$log.scale && nsc > 1) scale <- exp(scale)
          
              if(!is.element(model, c("Gumbel","Exponential"))) {

                X.sh <- designs$X.sh
	        nsh <- ncol(X.sh)
                shape <- rowSums(matrix(p[(nloc+nsc+1):np], x$n, nsh, byrow=TRUE) * X.sh)

              } else shape <- NULL

            } else {

              # blocks$designs should already exist
              if(is.element(model, c("PP","GEV","Gumbel","Weibull","Frechet"))) {
                X.loc <- x$blocks$designs$X.loc
	        nloc <- ncol(X.loc)
                loc <- rowSums(matrix(p[1:nloc], x$blocks$nBlocks, nloc, byrow=TRUE) * X.loc)

              } else loc <- NULL
          
              X.sc <- x$blocks$designs$X.sc
              nsc <- ncol(X.sc)
              scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], x$blocks$nBlocks, nsc, byrow=TRUE) * X.sc)
              if(x$par.models$log.scale) scale <- exp(scale)
          
              if(!is.element(model, c("Gumbel","Exponential"))) {

                X.sh <- x$blocks$designs$X.sh
	        nsh <- ncol(X.sh)
                shape <- rowSums(matrix(p[(nloc+nsc+1):np], x$blocks$nBlocks, nsh, byrow=TRUE) * X.sh)

              } else shape <- NULL

            }
	} else {

	    nr <- nrow(qcov)

	    if(is.element(model, c("GP", "Beta", "Exponential", "Pareto"))) nloc <- 0
	    else if(is.element("location", pnames)) nloc <- 1
	    else nloc <- sum(substring(pnames, 1, 2) == "mu")

	    if(x$par.models$log.scale) {

		if(!any(substring(pnames, 1, 3) == "phi")) nsc <- 1
		else nsc <- sum(substring(pnames, 1, 3) == "phi")

	    } else {

	        if(is.element("scale", pnames)) nsc <- 1
		else nsc <- sum(substring(pnames, 1, 3) == "sig")

	    }

	    if(is.element("shape", pnames)) nsh <- 1
	    else nsh <- sum(substring(pnames, 1, 2) == "xi")

	    if(!is.element(model, c("GP", "Beta", "Pareto"))) loc <- rowSums(matrix(rep(p[1:nloc], nr), nloc, nr, byrow = TRUE) * qcov[,1:nloc])
	    else loc <- NULL

	    scale <- rowSums(matrix(rep(p[(nloc + 1):(nloc + nsc)], nr), nsc, nr, byrow=TRUE) * qcov[,(nloc + 1):(nloc + nsc)])
	    if(x$par.models$log.scale && nsc > 1) scale <- exp(scale)

	    if(!is.element(model, c("Gumbel", "Exponential"))) shape <- rowSums(matrix(rep(p[(nloc + nsc + 1):(nloc + nsc + nsh)], nr), nsh, nr, byrow = TRUE) *
		qcov[,(nloc + nsc + 1):(nloc + nsc + nsh)])
	    else shape <- NULL
	    
	} # end of if else 'qcov' is null stmts.

    } # end of if '!tform' stmts.

    return(list(location=loc, scale=scale, shape=shape))
} # end of 'findpars.fevd.bayesian' function.

findpars.fevd.mle <- function(x, use.blocks=FALSE, ..., qcov = NULL) { # Eric 8/8/13 -- added 'qcov' argument.

    type <- tolower(x$type)

    p <- x$results$par
    np <- length(p)

    if(is.null(qcov)) { # Eric 8/8/13

        if(is.null(x$blocks) || !use.blocks) {  # CJP2

          n <- x$n  
          designs <- setup.design(x)

        } else {

          n <- x$blocks$nBlocks
          designs <- x$blocks$designs

        }
    
        if(is.element(type, c("gev", "weibull", "gumbel", "frechet","pp"))) {

            X.loc <- designs$X.loc
            nloc <- ncol(X.loc)
            loc <- rowSums(matrix(p[1:nloc], n, nloc, byrow=TRUE)*X.loc)

        } else {

            nloc <- 0
            loc <- NULL

        }

        X.sc <- designs$X.sc
        nsc <- ncol(X.sc)
        scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE)*X.sc)
        if(x$par.models$log.scale) scale <- exp(scale)

        if(!is.element(type, c("gumbel","exponential"))) {

            X.sh <- designs$X.sh
            nsh <- ncol(X.sh)
            shape <- rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], n, nsh, byrow=TRUE)*X.sh)

        } else {

            nsh <- 0
            shape <- 0

        }
    } else {
	
	pnames <- names(p)
	model <- x$type

	nr <- nrow(qcov)

        if(is.element(model, c("GP", "Beta", "Exponential", "Pareto"))) nloc <- 0
        else if(is.element("location", pnames)) nloc <- 1
        else nloc <- sum(substring(pnames, 1, 2) == "mu")

        if(x$par.models$log.scale) {

        if(!any(substring(pnames, 1, 3) == "phi")) nsc <- 1
        else nsc <- sum(substring(pnames, 1, 3) == "phi")

        } else {

            if(is.element("scale", pnames)) nsc <- 1
            else nsc <- sum(substring(pnames, 1, 3) == "sig")

        }

        if(is.element("shape", pnames)) nsh <- 1
        else nsh <- sum(substring(pnames, 1, 2) == "xi")

        if(!is.element(model, c("GP", "Beta", "Pareto"))) loc <- rowSums(matrix(rep(p[1:nloc], nr), nloc, nr, byrow = TRUE) * qcov[,1:nloc])
        else loc <- NULL

        scale <- rowSums(matrix(rep(p[(nloc + 1):(nloc + nsc)], nr), nsc, nr, byrow=TRUE) * qcov[,(nloc + 1):(nloc + nsc)])
        if(x$par.models$log.scale) scale <- exp(scale)

        if(!is.element(model, c("Gumbel", "Exponential"))) shape <- rowSums(matrix(rep(p[(nloc + nsc + 1):(nloc + nsc + nsh)], nr), nsh, nr, byrow = TRUE) *
                qcov[,(nloc + nsc + 1):(nloc + nsc + nsh)])
        else shape <- NULL

    } # end of if else 'qcov' is null stmts.

    return(list(location=loc, scale=scale, shape=shape))

} # end of 'findpars.fevd.mle' function.

profliker <- function(object, type=c("return.level", "parameter"), xrange=NULL, return.period=100, which.par=1, nint=20, plot=TRUE,
			    gr=NULL, method="BFGS", lower=-Inf, upper=Inf, control=list(), ...) {


    type <- tolower(type)
    type <- match.arg(type)
    designs <- setup.design(object)
    tmp <- datagrabber(object)

    if(object$data.name[2] == "") {
	y <- tmp
	data <- NULL
    } else {
	y <- tmp[,1]
	data <- tmp[,-1]
    }

    thresh <- object$threshold
    pars <- object$results$par
    pnames <- names(pars)

    args <- list(...)

    if(missing(nint) && object$type=="PP") nint <- 20 

    if(object$type=="PP" && type=="return.level") {
	object$type <- "GEV"
	plot.msg <- "Return levels based on equivalent GEV df"
    } else if(plot) plot.msg <- NULL

    if(type=="return.level") {

	if(!is.fixedfevd(object)) stop("profliker: Sorry, currently no support for profile likelihood with effective return levels.")
	if(is.null(xrange)) {

	    rl.hat <- rl.fevd(object, period=return.period)
	    if(rl.hat != 0) xrange <- c(rl.hat - 2 * abs(log2(abs(rl.hat))), rl.hat + 2 * abs(log2(abs(rl.hat))))
	    else xrange <- range(y, finite=TRUE)

	}

    } else {

        fv <- pars[which.par]
	if(is.null(xrange)) {

	    tmp <- summary(object, silent=TRUE)
	    tmp <- rbind(tmp$par, tmp$se.theta)

	    if(is.matrix(tmp)) {

		look <- c(tmp[,which.par])
		xrange <- c(look[1] - 4 * look[2], look[1] + 4 * look[2])
		th.names <- colnames(tmp)

	    } else {

		look <- tmp[which.par]
		xrange <- c(look - 4 * abs(log2(abs(look))), look + 4 * abs(log2(abs(look))))
		th.names <- names(tmp)

	    }

	    if((th.names[which.par] == "scale") && !object$par.models$log.scale) {

		xrange[1] <- max(xrange[1], 1e-10)
		if(xrange[1] > xrange[2]) xrange[2] <- 2 * xrange[1]

	    }

	} # end of if find 'xrange' stmts.

    } # end of if else 'return.levels' stmts.

    fixedvals <- seq(xrange[1], xrange[2],, nint)
    res <- numeric(nint) + NA

    if(type=="parameter") ipars <- pars[-which.par]
    else {

	if(is.element(object$type, c("GEV", "Gumbel", "Weibull", "Frechet"))) {

	    if(is.element("location", pnames)) ipars <- pars[pnames != "location"]
	    else {

		id <- substring(pnames, 1, 2) == "mu"
		ipars <- pars[!id]

	    }

	} else {

	    if(is.element("scale", pnames)) ipars <- pars[pnames != "scale"]
	    else if(is.element("log.scale", pnames)) ipars <- pars[pnames != "log.scale"]
	    else {

		istr <- substring(pnames, 1, 3)
		id <- (istr == "sig") | (istr == "phi")
		ipars <- pars[!id]

	    }

	}

    }

    for(i in 1:nint) {

	hold <- optim(par = ipars, fn = oevd.profpar, o=object, des=designs, x=y, data=data, u=thresh,
				fixed.index=which.par, fixed.value=fixedvals[i], which.type=type, rperiod=return.period,
				gr=gr, method=method, lower=lower, upper=upper, control=control)

 	ipars <- hold$par
	res[ i ] <- hold$value

    } # end of for 'i' loop.
   
     res <- -res 

     if(plot) {

	if(type == "parameter") xl <- names(pars)[which.par]
	else xl <- paste(return.period, "-", object$period.basis, " return level", sep="")

	xl <- paste(xl, plot.msg, sep="\n")

	if(is.null(args$main)) {

	    msg <- deparse(object$call)

	    if(is.null(args$xlab)) plot(fixedvals, res, type="l", xlab=xl, ylab="Profile Log-Likelihood", main=msg, ...)
	    else plot(fixedvals, res, type="l", ylab="Profile Log-Likelihood", main=msg, ...)

	} else {

	    if(is.null(args$xlab)) plot(fixedvals, res, type="l", xlab=xl, ylab="Profile Log-Likelihood", ...)
	    else plot(fixedvals, res, type="l", ylab="Profile Log-Likelihood", ...)

	}

	if(type=="parameter") abline(v=pars[which.par], lty=2, lwd=1.5)
	else abline(v=rl.fevd(object, period=return.period), lty=2, lwd=1.5)

    } # end of if 'plot' stmts.

    invisible(res)

} # end of 'profliker' function.

oevd.profpar <- function(p, o, des, x, data=NULL, u=NULL, fixed.index, fixed.value, which.type=c("return.level","parameter"), rperiod=100) {
    ##
    ## Function to be optimized by 'optim'.
    ## Takes the parameter design matrices and
    ## calculates the resulting vector of parameters,
    ## then calls 'levd' to get the negative log-likelihood.
    ## The only difference between this function and that of
    ## 'oevd' is that one of the parameters is held fixed in
    ## order to minimize a profile negative log-likelihood.
    ## Called by 'profliker' when a parameter's profile likelihood
    ## is of interest.
    which.type <- match.arg(which.type)
    type <- o$type
    span <- o$span
    npy <- o$npy
    n <- length(x)
    phi <- o$par.models$log.scale

    np <- length(p)+1

    if(which.type == "parameter") {

        p.tmp <- numeric(np) + NA
        p.tmp[-fixed.index] <- p
        p.tmp[fixed.index] <- fixed.value

    } else p.tmp <- p

    if(which.type=="parameter") {

        if(!is.element(type, c("GP","Beta","Pareto","Exponential"))) {

	    if(which.type=="parameter") {

                X.loc <- des$X.loc
                nloc <- ncol(X.loc)
                loc <- rowSums(matrix(p.tmp[1:nloc], n, nloc, byrow=TRUE)*X.loc)

	    } else nloc <- 1

        } else {

            nloc <- 0
            loc <- 0

        }

    } else {

	nloc <- 0
	loc <- 0

    } # end of if 'which.type' parameter stmts.

    if(which.type=="parameter" || is.element(type, c("GEV", "Gumbel", "Weibull", "Frechet", "PP"))) {

        X.sc <- des$X.sc
        nsc <- ncol(X.sc)
        if(phi) scale <- exp(rowSums(matrix(p.tmp[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE)*X.sc))
        else scale <- rowSums(matrix(p.tmp[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE)*X.sc)

    } else nsc <- 0

    if(!is.element(type, c("Gumbel","Exponential"))) {

        X.sh <- des$X.sh
        nsh <- ncol(X.sh)
        shape <- rowSums(matrix(p.tmp[(nloc+nsc+1):(nloc+nsc+nsh)], n, nsh, byrow=TRUE)*X.sh)

    } else {

        nsh <- 0
        shape <- 0

    }

    if(!is.element(type, c("GP","Beta","Pareto","Exponential")) && (which.type == "return.level")) {

	if(type != "Gumbel") loc <- fixed.value + (scale/shape) * (1 - (-log(1 - 1/rperiod))^(-shape))
	else loc <- fixed.value + scale * log(-log(1 - 1/rperiod))

    } else if(which.type == "return.level") {

	lam <- mean(x > u)
	m <- rperiod * npy
	if(type != "Exponential") {
	    zid <- shape == 0
	    if(any(zid)) shape[zid] <- 1e-10
	    scale <- ((fixed.value - u) * shape)/((m * lam)^shape - 1)

	} else scale <- (fixed.value - u)/log(m * lam)
    }

    if(is.null(span)) res <- levd(x=x, threshold=u, location=loc, scale=scale, shape=shape, type=type, npy=npy, infval=1e10)  
    else res <- levd(x=x, threshold=u, location=loc, scale=scale, shape=shape, type=type, span=span, npy=npy, infval=1e10)  
    # if(class(res) == "try-error") res <- 1e10

    return(res)

} # end of 'oevd.profpar' function.

rl.fevd <- function(x, period=100) {

    pmods <- x$par.models

    if(is.fixedfevd(x)) {

	p <- x$results$par

	if(any(names(p) == "log.scale")) {

	    id <- names(p) == "log.scale"
	    p[id] <- exp(p[id])
	    names(p)[id] <- "scale"

	}

	if(!is.null(x$threshold)) lam <- mean(datagrabber(x)[,1] > x$threshold)
	else lam <- NULL

	return(rlevd(period=period, loc=p["location"], scale=p["scale"], shape=p["shape"], threshold=x$threshold, type=x$type, npy=x$npy, rate=lam))

    } else return(erlevd(x, period=period))

} # end of 'rl.fevd' function.

rextRemes <- function(x, n, ...) {

    UseMethod("rextRemes", x)

} # end of 'rextRemes' function.

rextRemes.fevd <- function(x, n, ...) {

    if(x$method=="GMLE") newcl <- "mle"
    else newcl <- tolower(x$method)
    class(x) <- paste("fevd.", newcl, sep="")
    UseMethod("rextRemes", x)

} # end of 'rextRemes.fevd' function. 

rextRemes.fevd.lmoments <- function(x, n, ...) {

    p <- x$results
    pnames <- names(p)

    if(any(pnames == "location")) loc <- p["location"]
    else loc <- 0

    scale <- p["scale"]

    if(any(pnames == "shape")) shape <- p["shape"]
    else shape <- 0

    type <- x$type

    if(is.element(type, c("PP","GP","Exponential","Beta","Pareto"))) u <- x$threshold
    else u <- 0
 
    if(type=="PP") {
	scale <- scale + shape * (u - loc)
	type <- "GP"
    }

    if(is.element(type, c("Gumbel","Weibull","Frechet"))) type <- "GEV"
    else if(is.element(type, c("Beta","Exponential","Pareto"))) type <- "GP"

    out <- revd(n=n, loc=loc, scale=scale, shape=shape, threshold=u, type=type)

    return(out)

} # end of 'rextRemes.fevd.lmoments' function.

rextRemes.fevd.bayesian <- function(x, n, ..., burn.in=499, FUN="mean", qcov = NULL) {

    type <- x$type

    f <- match.fun(FUN)
    
    if(is.fixedfevd(x)) {
 
    p <- x$results
    np <- ncol(p) - 1
    p <- p[,1:np]
    pnames <- colnames(p)

    if(burn.in > 0) p <- p[-(1:burn.in),]
      
    if(is.element("log.scale",pnames)) {

            id <- pnames == "log.scale"
            p[,id] <- exp(p[,id])
            pnames[id] <- "scale"
            colnames(p) <- pnames

    }

    if(FUN == "mean") p <- colMeans(p)
    else p <- apply(p, 2, f)

    if(is.element(type, c("PP","GEV","Gumbel","Weibull","Frechet"))) {

        loc <- p["location"]
        nloc <- 1

    } else loc <- 0

    scale <- p["scale"]

    if(!is.element(type, c("Gumbel","Exponential"))) shape <- p["shape"]
    else shape <- 0

    if(type == "PP") {

       scale <- scale + shape * (x$threshold - loc)
       type <- "GP"

    }

    return(revd(n=n, loc=loc, scale=scale, shape=shape, threshold=x$threshold, type=type))

    } else {

        # ytrans <- trans(x)
        if(missing(n)) n <- 1
        p <- findpars(x, burn.in=burn.in, FUN=FUN, qcov = qcov)
    
        if(is.null(qcov)) K <- x$n
        else K <- dim(qcov)[1]

        if(!is.null(qcov)) u <- qcov[,"threshold"]  
        else if(is.null(x$threshold)) u <- rep(0, K)
        else u <- x$threshold

        if(is.null(p$location)) loc <- rep(0, K)
        else loc <- p$location
      
        if(is.null(p$shape)) {

          if(type=="Gumbel") shape <- rep(0, K)
          else shape <- rep(1e-10, K)

      } else shape <- p$shape
      
      if(type != "PP") theta <- cbind(u, loc, p$scale, shape)
      else theta <- cbind(u, loc, p$scale + shape * (u - loc), shape)

      if(is.null(qcov)) y <- c(datagrabber(x, cov.data=FALSE))

      if(is.element(type,c("GEV","Gumbel","Weibull","Frechet"))) {

        out <- matrix(NA, K, n)
        if(is.null(qcov)) o <- order(y)

        for(i in 1:n) {

          z <- revd(K, type="GEV")
          if(is.null(qcov)) z <- sort(z)[o]

          out[,i] <- revtrans.evd(z=z, threshold=theta[,1], location=theta[,2], scale=theta[,3], shape=theta[,4], type=type)

        } # end of for 'i' loop.

      } else {

	if(is.null(qcov)) {

            eid <- y > u
            m <- sum(eid)
	    theta <- theta[eid,]
            # y <- y[eid]
            # o <- order(y) # Eric 8/8/13 -- looks like this was not used.

	} else {

	    m <- K
	    
	}

	out <- matrix(NA, m, n)

        for(i in 1:n) {

          z <- revd(m, type="GP")
          out[,i] <- revtrans.evd(z=z, threshold=theta[,1], location=theta[,2], scale=theta[,3], shape=theta[,4], type="GP")

        } # end of for 'i' loop.

      } # end of if else type of EVD stmts.

      return(out)

    } # end of if 'is.fixedfevd' stmts.

} # end of 'rextRemes.fevd.bayesian' function.

rextRemes.fevd.mle <- function(x, n, ..., qcov = NULL) {

  type <- x$type
  
  if(is.fixedfevd(x)) {

    if(missing(n)) n <- x$n
    p <- x$results$p
    pnames <- names(p)
    if(is.element("log.scale", pnames)) {
      id <- pnames == "log.scale"
      p[id] <- exp(p[id])
      pnames[id] <- "scale"
      names(p) <- pnames
    }
    scale <- p["scale"]
    if(x$par.models$log.scale) scale <- exp(scale)
    if(type=="GEV") {
      loc <- p["location"]
      shape <- p["shape"]
      u <- 0
    } else if(type=="GP") {
      loc <- 0
      u <- x$threshold
      shape <- p["shape"]
    } else if(type=="PP") {
      loc <- p["location"]
      u <- x$threshold
      shape <- p["shape"]
      scale <- scale + shape * (u - loc)
      type <- "GP"
    } else if(type=="Gumbel") {
      loc <- p["location"]
      shape <- 0
      type <- "GEV"
      u <- 0
    } else if(is.element(type, c("Weibull","Frechet"))) {
      loc <- p["location"]
      shape <- p["shape"]
      u <- 0
      type <- "GEV"
    } else if(type=="Exponential") {
      loc <- 0
      shape <- 0
      u <- x$threshold
      type <- "GP"
    } else if(is.element(type, c("Beta","Pareto"))) {
      loc <- 0
      shape <- p["shape"]
      u <- x$threshold
      type <- "GP"
    }
    
    return(revd(n=n, loc=loc, scale=scale, shape=shape, threshold=u, type=type))
    
  } else {

    if(missing(n)) n <- 1
    p <- findpars(x, qcov = qcov)
   
    if(is.null(qcov)) K <- x$n
    else K <- dim(qcov)[1]

    if(!is.null(qcov)) u <- qcov[,"threshold"] 
    else if(is.null(x$threshold)) u <- rep(0, K)
    else u <- x$threshold
    
    if(is.null(p$location)) loc <- rep(0, K)
    else loc <- p$location
    
    if(is.null(p$shape)) {

      if(type=="Gumbel") shape <- rep(0, K)
      else shape <- rep(1e-10, K)

    } else shape <- p$shape
    
    scale <- p$scale
    if(type=="PP") scale <- scale + shape * (u - loc)
    
    theta <- cbind(u, loc, scale, shape)

    if(is.null(qcov)) y <- c(datagrabber(x, cov.data=FALSE))

    if(is.element(type,c("GEV","Gumbel","Weibull","Frechet"))) {

      out <- matrix(NA, K, n)

      for(i in 1:n) {

        z <- revd(K, type="GEV")
        out[,i] <- revtrans.evd(z=z, threshold=theta[,1], location=theta[,2], scale=theta[,3], shape=theta[,4], type=type)

      } # end of for 'i' loop.

    } else {

      if(is.null(qcov)) {

	  eid <- y > u
          m <- sum(eid)
	  theta <- theta[eid,]

      } else {

	  m <- K
      }

      out <- matrix(NA, m, n)

      if(type=="PP") type2 <- "GP" # Eric -- 8/8/13 -- should we change this to GEV?
      else type2 <- type

      for(i in 1:n) {

        z <- revd(m, type="GP")
        out[,i] <- revtrans.evd(z=z, threshold=theta[,1], location=theta[,2], scale=theta[,3], shape=theta[,4], type=type2)

      } # end of for 'i' loop.
    }
    return(out)
  }
} # end of 'rextRemes.fevd' function.

pextRemes <- function(x, q, lower.tail=TRUE, ...) {
    if(missing(q)) stop("pextRemes: Must specify q argument.")
    UseMethod("pextRemes", x)
} # end of 'pextRemes' function.

pextRemes.fevd <- function(x, q, lower.tail=TRUE, ..., qcov=NULL) {
    if(x$method=="GMLE") newcl <- "mle"
    else newcl <- tolower(x$method)
    class(x) <- paste("fevd.", newcl, sep="")
    UseMethod("pextRemes", x)
} # end of 'pextRemes.fevd' function.

pextRemes.fevd.lmoments <- function(x, q, lower.tail=TRUE, ...) {

    type <- x$type

    p <- x$results
    pnames <- names(p)

    if(is.element("location", pnames)) loc <- p["location"]
    else loc <- 0

    scale <- p["scale"]

    if(is.element("shape", pnames)) shape <- p["shape"]
    else shape <- 0

    if(!is.null(x$threshold)) u <- x$threshold
    else u <- 0

    if(type=="Gumbel") type <- "GEV"
    else if(type=="Exponential") type <- "GP"
    else if(type=="PP") {
        scale <- scale + shape * (u - loc)
        type <- "GP"
    } else if(is.element(type, c("Weibull","Frechet"))) type <- "GEV"
    else if(is.element(type, c("Beta","Pareto"))) type <- "GP"

    return(pevd(q=q, loc=loc, scale=scale, shape=shape, threshold=u, lambda=x$rate, npy=x$npy, type=type, lower.tail=lower.tail))

} # end of 'pextRemes.fevd.lmoments' function.

pextRemes.fevd.bayesian <- function(x, q, lower.tail=TRUE, ..., qcov=NULL, burn.in=499, FUN="mean") {

    type <- x$type

    p <- x$results
    np <- ncol(p) - 1
    p <- p[,1:np]
    pnames <- colnames(p)
    if(burn.in > 0) p <- p[-(1:burn.in),]

    f <- match.fun(FUN)
    p <- apply(p, 2, f)

    phi <- x$par.models$log.scale

    if(is.fixedfevd(x)) {

	if(is.element("location", pnames)) loc <- p["location"]
	else loc <- 0

	if(phi) scale <- exp(p["log.scale"])
	else scale <- p["scale"]

	if(is.element("shape", pnames)) shape <- p["shape"]
	else shape <- 0

	if(!is.null(x$threshold)) u <- x$threshold
	else u <- 0

	if(type=="Gumbel") type <- "GEV"
	else if(type=="Exponential") type <- "GP" 
	else if(type=="PP") {
	    scale <- scale + shape * (u - loc)
	    type <- "GP"
	} else if(is.element(type, c("Weibull","Frechet"))) type <- "GEV"
	else if(is.element(type, c("Beta","Pareto"))) type <- "GP"

	return(pevd(q=q, loc=loc, scale=scale, shape=shape, threshold=u, lambda=x$rate, npy=x$npy, type=type, lower.tail=lower.tail))

    } else {

	n <- length(q)

	if(is.element("mu", substring(pnames, 1, 2))) nloc <- sum(substring(pnames, 1, 2) == "mu")
        else if(is.element("location", pnames)) nloc <- 1
        else nloc <- 0

        if(is.element("phi", substring(pnames, 1, 3))) nsc <- sum(substring(pnames, 1, 3) == "phi")
        else if(is.element("sig", substring(pnames, 1, 3))) nsc <- sum(substring(pnames, 1, 3) == "sig")
        else nsc <- 1

        if(is.element("shape", pnames)) nsh <- 1
        else if(is.element("xi", substring(pnames, 1, 2))) nsh <- sum(substring(pnames, 1, 2) == "xi")
        else nsh <- 0

	if(is.null(qcov)) {
            warning("pextRemes: qcov argument not supplied to non-stationary model.  Using intercept terms only and/or if non-constant threshold, using first value.")
            qcov <- matrix(0, n, np+1)
            qcov[,1] <- 1
            if(nloc > 0) qcov[,nloc+1] <- 1
            qcov[,nloc+nsc+1] <- 1
            if(!is.null(x$threshold)) qcov[,np+1] <- x$threshold[1]
	    else qcov[,np + 1] <- 0
            colnames(qcov) <- c(pnames, "threshold")
        }

        if(!is.matrix(qcov)) stop("pextRemes: invalid type for qcov argument.  Must be a matrix.")

	if(dim(qcov)[2] == np+1) {
            if(is.null(colnames(qcov))) colnames(qcov) <- c(pnames, "threshold")
            else if(!is.element("threshold", colnames(qcov))) {
                warning("pextRemes: qcov columns are named and dimension suggests threshold, but no threshold name.  Assuming last column is threshold.")
                colnames(qcov) <- c(pnames, "threshold")
            }
        } # end of if 'dim(qcov)[2] == np + 1' stmts.

	u <- qcov[,"threshold"]

        if(nloc==0) loc <- rep(0, n)
        else loc <- rowSums(matrix(p[1:nloc], n, nloc, byrow=TRUE) * qcov[,1:nloc])

        if(nsh == 0) {
            if(type=="Gumbel") shape <- rep(0, n)
            else shape <- rep(1e-10, n)
        } else shape <- rowSums(matrix(p[(nloc+nsc+1):np], n, nsh, byrow=TRUE) * qcov[,(nloc+nsc+1):np])

        scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE) * qcov[,(nloc+1):(nloc+nsc)])
        if(phi) scale <- exp(scale)

	if(type=="PP") {
	    scale <- scale + shape * (u - loc)
	    type <- "GP"
	} else if(is.element(type, c("Gumbel","Weibull","Frechet"))) type <- "GEV"
	else if(is.element(type, c("Exponential","Beta","Pareto"))) type <- "GP"

	theta <- cbind(q, u, loc, scale, shape)

        pfun <- function(th, type, rate, npy, lower.tail){
	    return(pevd(q=th[1], loc=th[3], scale=th[4], shape=th[5], threshold=th[2], lambda=rate, npy=npy, lower.tail=lower.tail))
	} # end of internal 'pfun' function.

        out <- apply(theta, 1, pfun, type=type, rate=x$rate, npy=x$npy, lower.tail=lower.tail)
        
        return(out)

    } # end of if else 'stationary model' stmts.

} # end of 'pextRemes.fevd.bayesian' function.

pextRemes.fevd.mle <- function(x, q, lower.tail=TRUE, ..., qcov=NULL) {

    type <- x$type

    p <- x$results$par
    pnames <- names(p)
    np <- length(p)
    phi <- x$par.models$log.scale

    if(is.fixedfevd(x)) {

        if(phi) scale <- exp(p["log.scale"])
        else scale <- p["scale"]

        if(type=="GEV") {
            loc <- p["location"]
            shape <- p["shape"]
            u <- 0
        } else if(type=="GP") {
            loc <- 0
            u <- x$threshold
            shape <- p["shape"]
        } else if(type=="PP") {
            loc <- p["location"]
            u <- x$threshold
            shape <- p["shape"]
            scale <- scale + shape * (u - loc)
            type <- "GP"
        } else if(type=="Gumbel") {
            loc <- p["location"]
            shape <- 0
            type <- "GEV"
            u <- 0
        } else if(is.element(type, c("Weibull","Frechet"))) {
            loc <- p["location"]
            shape <- p["shape"]
            u <- 0
            type <- "GEV"
        } else if(type=="Exponential") {
            loc <- 0
            shape <- 0
            u <- x$threshold
            type <- "GP"
        } else if(is.element(type, c("Beta","Pareto"))) {
            loc <- 0
            shape <- p["shape"]
            u <- x$threshold
            type <- "GP"
        }

        return(pevd(q=q, loc=loc, scale=scale, shape=shape, threshold=u, lambda=x$rate, npy=x$npy, type=type, lower.tail=lower.tail))

    }  else {

	n <- length(q)

	if(is.element("mu", substring(pnames, 1, 2))) nloc <- sum(substring(pnames, 1, 2) == "mu")
        else if(is.element("location", pnames)) nloc <- 1
        else nloc <- 0

	if(is.element("phi", substring(pnames, 1, 3))) nsc <- sum(substring(pnames, 1, 3) == "phi")
	else if(is.element("sig", substring(pnames, 1, 3))) nsc <- sum(substring(pnames, 1, 3) == "sig")
	else nsc <- 1

	if(is.element("shape", pnames)) nsh <- 1
	else if(is.element("xi", substring(pnames, 1, 2))) nsh <- sum(substring(pnames, 1, 2) == "xi")
	else nsh <- 0

	if(is.null(qcov)) {
	    warning("pextRemes: qcov argument not supplied to non-stationary model.  Using intercept terms only and/or if non-constant threshold, using first value.")
	    qcov <- matrix(0, n, np+1)
	    qcov[,1] <- 1
	    if(nloc > 0) qcov[,nloc+1] <- 1
	    qcov[,nloc+nsc+1] <- 1
	    if(!is.null(x$threshold)) qcov[,np+1] <- x$threshold[1]
	    else qcov[,np + 1] <- 0
	    colnames(qcov) <- c(pnames, "threshold")
	}

	if(!is.matrix(qcov)) stop("pextRemes: invalid type for qcov argument.  Must be a matrix.")

	if(dim(qcov)[2] == np+1) {
	    if(is.null(colnames(qcov))) colnames(qcov) <- c(pnames, "threshold")
	    else if(!is.element("threshold", colnames(qcov))) {
		warning("pextRemes: qcov columns are named and dimension suggests threshold, but no threshold name.  Assuming last column is threshold.")
		colnames(qcov) <- c(pnames, "threshold")
	    }
	}

	if(!is.element("threshold", colnames(qcov))) {
	    if(is.null(x$threshold)) qcov <- cbind(qcov, 0)
	    else if(length(x$threshold)==1) qcov <- cbind(qcov, x$threshold)
	    else {
		warning("pextRemes: non-constant threshold, but threshold values not supplied to qcov argument.  Using first value only")
		qcov <- cbind(qcov, x$threshold[1])
	    }
	    colnames(qcov) <- c(pnames, "threshold")
	}

	if((dim(qcov)[1] != n) || (!is.element(dim(qcov)[2], c(np,np+1)))) {
	    stop("pextRemes: invalid qcov argument.  Must be a matrix with rows equal to the length of q and columns equal to the number of parameters, np, or np + 1.")
	}

	u <- qcov[,"threshold"]

        if(nloc==0) loc <- rep(0, n)
        else loc <- rowSums(matrix(p[1:nloc], n, nloc, byrow=TRUE) * qcov[,1:nloc])

        if(nsh == 0) {
            if(type=="Gumbel") shape <- rep(0, n)
            else shape <- rep(1e-10, n)
        } else shape <- rowSums(matrix(p[(nloc+nsc+1):np], n, nsh, byrow=TRUE) * qcov[,(nloc+nsc+1):np])

	scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE) * qcov[,(nloc+1):(nloc+nsc)])
	if(phi) scale <- exp(scale)

        if(type=="PP") {
	    scale <- scale + shape * (u - loc)
	    type <- "GP"
	} else if(is.element(type, c("Gumbel","Weibull","Frechet"))) type <- "GEV"
	else if(is.element(type, c("Exponential","Beta","Pareto"))) type <- "GP"

        theta <- cbind(q, u, loc, scale, shape)
	pfun <- function(th, type, rate, npy, lower.tail) {
	    return(pevd(q=th[1], loc=th[3], scale=th[4], shape=th[5], threshold=th[2], lambda=rate, npy=npy, type=type, lower.tail=lower.tail))
	} # end of internal 'pfun' function.
	out <- apply(theta, 1, pfun, type=type, rate=x$rate, npy=x$npy, lower.tail=lower.tail)
        
        return(out)

    } # end of if else stationary model stmts.

} # end of 'pextRemes.fevd.mle' function

is.fixedfevd <- function(x) {

    p <- x$par.models
    return(all(c(check.constant(x$threshold), check.constant(p$threshold), check.constant(p$location), check.constant(p$scale), check.constant(p$shape))))

} # end of 'is.fixedfevd' function.

ci.fevd <- function(x, alpha=0.05, type=c("return.level", "parameter"), return.period=100, which.par, R=502, ...) {

    if(x$method == "GMLE") newcl <- "fevd.mle"
    else newcl <- paste("fevd.", tolower(x$method), sep="")
    class(x) <- newcl
    UseMethod("ci",x)

} # end of 'ci.fevd' function.

ci.fevd.bayesian <- function(x, alpha=0.05, type=c("return.level", "parameter"), return.period=100, which.par=1, FUN="mean", burn.in=499, tscale=FALSE, ...) {

    type <- tolower(type)
    type <- match.arg(type)

    f <- match.fun(FUN)

    const <- is.fixedfevd(x)
    if(type=="return.level" && !const) return(ci.rl.ns.fevd.bayesian(x = x, alpha = alpha, return.period = return.period, FUN = FUN, burn.in = burn.in, ...))

    pars <- x$results
    np <- dim(pars)[2] - 1
    pars <- pars[,1:np]
    if(burn.in > nrow(pars)) stop("ci: burn.in > number of MCMC iterations (iter).  Must be < iter.")
    if(burn.in != 0) pars <- pars[-(1:burn.in),]
    pnames <- colnames(pars)
    iters <- dim(pars)[1]

    if(type=="parameter" && missing(which.par)) which.par <- 1:np

    if(is.element("log.scale",pnames)) {
	id <- pnames=="log.scale"
	pars[,id] <- exp(pars[,id])
	pnames[id] <- "scale"
	colnames(pars) <- pnames
    }

    conf.level <- paste(round((1 - alpha)*100, digits=2), "%", sep="")

    if(type == "return.level") {

	if(is.element("location",pnames)) loc <- pars[,"location"]
	else loc <- rep(0, iters)

	scale <- pars[,"scale"]

	if(is.element("shape",pnames)) shape <- pars[,"shape"]
	else shape <- rep(0, iters)

	if(!is.null(x$threshold)) u <- x$threshold
	else u <- 0

	theta.sam <- cbind(loc, scale, shape)

	rlfun <- function(th, pd, u, type, npy, rate) rlevd(period=pd, loc=th[1], scale=th[2], shape=th[3], threshold=u, type=type, npy=npy, rate=rate)
	theta.sam <- apply(theta.sam, 1, rlfun, pd=return.period, u=u, type=x$type, npy=x$npy, rate=x$rate)

	if(is.matrix(theta.sam)) {

	    if(FUN=="mean") theta.hat <- rowMeans(theta.sam)
	    else theta.hat <- apply(theta.sam, 1, f, ...)
	    bds <- apply(theta.sam, 1, quantile, probs=c(alpha/2,1-alpha/2))
	    est.names <- paste(names(theta.hat), "-", x$period.basis, " level", sep="")
	    bds.names <- rownames(bds)
	    out <- cbind(bds[1,], theta.hat, bds[2,])
            rownames(out) <- est.names
            colnames(out) <- c(bds.names[1], "Posterior Mean", bds.names[2])

	} else {

	    if(FUN=="mean") theta.hat <- mean(theta.sam)
	    else theta.hat <- f(theta.sam)

	    bds <- quantile(theta.sam, probs=c(alpha/2,1-alpha/2))
	    out <- c(bds[1], theta.hat, bds[2])
	    names(out) <- c(paste(conf.level, " lower CI", sep=""),
			    paste("Posterior Mean ", return.period, "-", x$period.basis, " level", sep=""),
			    paste(conf.level, " upper CI", sep=""))
	}

    } else if(type=="parameter") {
	theta.hat <- colMeans(pars)
	est.names <- names(theta.hat)

	if(tscale) {

	    if(!const && !is.element("scale",est.names) && !is.element("shape",est.names)) stop("ci: invalid argument configurations.")
            if(!is.element(x$type, c("GP","Beta","Pareto"))) stop("ci: invalid argument configurations.")
            pars[,"scale"] <- pars[,"scale"] - pars[,"shape"] * x$threshold
            est.names[est.names == "scale"] <- "tscale"
            colnames(pars) <- est.names
            theta.hat["scale"] <- theta.hat["scale"] - theta.hat["shape"] * x$threshold
            names(theta.hat) <- est.names

        }

	theta.hat <- theta.hat[which.par]
	est.names <- est.names[which.par]

	bds <- apply(pars, 2, quantile, probs=c(alpha/2,1-alpha/2))
	if(is.matrix(bds)) bds <- bds[,which.par]
	else bds <- bds[which.par]

	if(is.matrix(bds)) {

	    bds.names <- rownames(bds)
	    out <- cbind(bds[1,], theta.hat, bds[2,])
	    rownames(out) <- est.names
	    colnames(out) <- c(bds.names[1], "Posterior Mean", bds.names[2])

	} else {

	    out <- c(bds[1], theta.hat, bds[2])
	    names(out) <- c(paste(conf.level, " lower CI", sep=""),
                            paste("Posterior Mean ", est.names, " parameter", sep=""),
			    paste(conf.level, " upper CI", sep=""))

	}

    } else stop("ci: invalid type parameter.")

    attr(out, "data.name") <- x$call
    attr(out, "method") <- "Quantiles of MCMC Sample from Posterior Distribution"
    attr(out, "conf.level") <- (1 - alpha) * 100
    class(out) <- "ci"
    return(out)
} # end of 'ci.fevd.bayesian' function.


ci.fevd.lmoments <- function(x, alpha=0.05, type=c("return.level", "parameter"), return.period=100, which.par, R=502, tscale=FALSE, return.samples=FALSE, ...) {

    type <- tolower(type)
    type <- match.arg(type)

    p <- x$results
    pnames <- names(p)
    np <- length(p)

    if(type=="parameter" && missing(which.par)) which.par <- 1:np

    if(is.element(x$type, c("GEV","Gumbel","Frechet","Weibull"))) n <- x$n
    else {
	look <- datagrabber(x, cov.data=FALSE)
        eid <- look > x$threshold
	n <- sum(eid)
    }

    z <- rextRemes(x, n=R * n)
    z <- matrix(z, n, R)

    if(!is.element(x$type, c("GEV","Gumbel","Frechet","Weibull"))) {
	z2 <- matrix(x$threshold, x$n, R)
	z2[eid,] <- z
	z <- z2
    }

    bfun <- function(y, ...) distill(fevd(y, ...))
    if(x$type=="GEV") sam <- apply(z, 2, bfun, method="Lmoments")
    else sam <- apply(z, 2, bfun, threshold=x$threshold, type=x$type, method="Lmoments", span=x$span)

    if(tscale) {
	if(!is.element(x$type, c("GP","Beta","Pareto"))) stop("ci: invalid argument configurations.")
	sam["scale",] <- sam["scale",] - sam["shape",] * x$threshold
	pnames[pnames == "scale"] <- "tscale"
	rownames(sam) <- pnames
    }

    if(return.samples) {
	if(is.matrix(sam)) out <- t(sam)
	else out <- cbind(sam)
	if(type=="parameter") return(out)
    }

    if(type=="return.level") {

	if(is.element("location",pnames)) {
	    loc <- sam["location",]
	    loc0 <- p["location"]
	} else loc <- loc0 <- 0

	scale <- sam["scale",]
	scale0 <- p["scale"]

	if(is.element("shape",pnames)) {
	    shape <- sam["shape",]
	    shape0 <- p["shape"]
	} else shape <- shape0 <- 0

	if(!is.null(x$threshold)) u <- x$threshold
	else u <- 0

	theta <- rbind(loc, scale, shape)
	rlfun <- function(th, ...) rlevd(period=return.period, loc=th[1], scale=th[2], shape=th[3], ...)
	
	sam <- apply(theta, 2, rlfun, threshold=u, type=x$type, npy=x$npy, rate=x$rate)

	if(is.matrix(sam)) rownames(sam) <- paste(return.period, "-", x$period.basis, sep="")
	else sammy.name <- paste(return.period, "-", x$period.basis, " level", sep="")
	
	if(return.samples) {
	    if(is.matrix(sam)) out <- cbind(out, t(sam))
	    else { 
		out <- cbind(out, c(sam))
		onames <- colnames(out)
		colnames(out) <- c(onames[-length(onames)], sammy.name) 
	    }
	    return(out)
	}

	theta.hat <- rlevd(period=return.period, loc=loc0, scale=scale0, shape=shape0, threshold=u, type=x$type, npy=x$npy, rate=x$rate)

    } else {
	theta.hat <- p
	if(tscale) {
            if(!is.element(x$type, c("GP","Beta","Pareto"))) stop("ci: scale parameter is not a function of threshold for this df.")
            sam["scale",] <- sam["scale",] - sam["shape",] * x$threshold
            pnames[pnames == "scale"] <- "tscale"
            rownames(sam) <- pnames
	    theta.hat["scale"] <- theta.hat["scale"] - theta.hat["shape"] * x$threshold
            names(theta.hat) <- pnames
        }

	sam <- sam[which.par,]
	if(!is.matrix(sam)) names(sam) <- pnames[which.par]
    }

    if(is.matrix(sam)) {
        out <- apply(sam, 1, quantile, probs=c(alpha/2, 1 - alpha/2))
        out.names <- rownames(out)
        out <- rbind(out[1,], theta.hat, out[2,])
        rownames(out) <- c(out.names[1], "Estimate", out.names[2])
        colnames(out) <- rownames(sam)
        out <- t(out)
    } else {
        out <- quantile(sam, probs=c(alpha/2, 1 - alpha/2))
	onames <- names(out)
        out <- c(out[1], theta.hat, out[2])
	if(type=="parameter") names(out) <- c(onames[1], pnames[which.par], onames[2])
	else names(out) <- c(onames[1], paste(return.period, "-", x$period.basis, " level", sep=""), onames[2])
    }
    attr(out, "method") <- "Parametric Bootstrap"
    attr(out, "data.name") <- x$call
    attr(out, "conf.level") <- (1 - alpha) * 100
    attr(out, "R") <- R
    class(out) <- "ci"
    return(out)
} # end of 'ci.fevd.lmoments' function.

ci.fevd.mle <- function(x, alpha=0.05, type=c("return.level", "parameter"), return.period=100, which.par=1, R=502,
		 	    method=c("normal","boot","proflik"), xrange=NULL, nint=20, verbose=FALSE, tscale=FALSE, 
			    return.samples=FALSE, ...) {

    if(missing(method)) miss.meth <- TRUE
    else miss.meth <- FALSE

    method <- tolower(method)
    method <- match.arg(method)

    type <- tolower(type)
    type <- match.arg(type)

    theta.hat <- x$results$par
    theta.names <- names(theta.hat)
    np <- length(theta.hat)

    if(type=="parameter" && missing(which.par)) which.par <- 1:np

    if(any(theta.names=="log.scale")) {
	id <- theta.names=="log.scale"
	theta.hat[id] <- exp(theta.hat[id])
	theta.names[id] <- "scale"
	names(theta.hat) <- theta.names
    }

    const <- is.fixedfevd(x)

    if(type=="return.level") par.name <- paste(return.period, "-", x$period.basis, " return level", sep="")
    else if(type=="parameter") par.name <- theta.names[which.par]

    if(type=="return.level" && !const) {
        return(ci.rl.ns.fevd.mle(x = x, alpha = alpha, return.period = return.period, method = method, verbose = verbose, return.samples = return.samples, ...))
# stop("ci: Sorry, no current functionality for finding CI of effective return levels.")
    }
 
    if(type=="parameter") p <- theta.hat[which.par]
    else {
	if(is.element(x$type, c("PP","GP","Beta","Pareto","Exponential"))) lam <- mean(c(datagrabber(x)[,1]) > x$threshold)
	else lam <- 1

	if(is.element(x$type, c("PP","GEV","Gumbel","Weibull","Frechet"))) loc <- theta.hat["location"]
        else loc <- 0
    	
	scale <- theta.hat["scale"]

        if(!is.element(x$type, c("Gumbel","Exponential"))) shape <- theta.hat["shape"]
        else shape <- 0

        if(x$type == "PP") mod <- "GEV"
        else mod <- x$type

	p <- rlevd(period=return.period, loc=loc, scale=scale, shape=shape, threshold=x$threshold, type=mod, npy=x$npy, rate=lam)
    }
 
    if(verbose) {
        cat("\n", "Preparing to calculate ", (1 - alpha)*100, "% CI for ",
                            ifelse(type=="return.level", paste(return.period, "-", x$period.basis, " return level", sep=""),
                                                         paste(par.name, " parameter", sep="")), "\n")
        cat("\n", "Model is ", ifelse(const, " fixed", " non-stationary."), "\n")
	if(method=="normal") cat("\n", "Using Normal Approximation Method.\n")
	else if(method=="boot") cat("\n", "Using Bootstrap Method.\n")
	else if(method=="proflik") cat("\n", "Using Profile Likelihood Method.\n")
    }

    if(method=="normal") {

	method.name <- "Normal Approx."
	z.alpha <- qnorm(alpha/2, lower.tail=FALSE)
	cov.theta <- parcov.fevd(x)
        if(is.null(cov.theta)) stop("ci: Sorry, unable to calculate the parameter covariance matrix.  Maybe try a different method.")
        var.theta <- diag(cov.theta)
	if(any(var.theta < 0)) stop("ci: negative Std. Err. estimates obtained.  Not trusting any of them.")

	if(type=="parameter") {

	    se.theta <- sqrt(var.theta)

	    if(tscale) {
                if(!const && !is.element("scale",theta.names) && !is.element("shape",theta.names) && !all(x$threshold == x$threshold[1])) {
		    stop("ci: invalid argument configurations.")
	 	}
                if(!is.element(x$type, c("GP","Beta","Pareto"))) stop("ci: invalid argument configurations.")
		theta.hat["scale"] <- theta.hat["scale"] - theta.hat["shape"] * x$threshold
                theta.names[theta.names == "scale"] <- "tscale"
		if(!any(theta.names[which.par] == "tscale")) stop("ci: invalid argument configurations.")
                names(theta.hat) <- theta.names
		p <- theta.hat[which.par]
		d <- rbind(1, -x$threshold)
		names(se.theta) <- theta.names
 		se.theta["tscale"] <- sqrt(t(d) %*% cov.theta %*% d)
            } else se.theta <- sqrt(var.theta)[which.par]

	    se.theta <- se.theta[which.par]
   	    par.name <- theta.names[which.par]
	} else if(type=="return.level") {

	    grads <- rlgrad.fevd(x, period=return.period)
	    grads <- t(grads)
	    if(is.element(x$type,c("GP","Beta","Pareto", "Exponential"))) {
		if(x$type=="Exponential") cov.theta <- diag(c(lam * (1 - lam)/x$n, var.theta))
		else cov.theta <- rbind(c(lam * (1 - lam)/x$n, 0, 0), cbind(0, cov.theta))
	    } else lam <- 1
	    
	    var.theta <- t(grads) %*% cov.theta %*% grads # CJP2: took out sqrt as taking sqrt of off-diags

	} else stop("ci: invalid type argument.  Must be return.level or parameter.")

	if(length(p) > 1) {
	    if(type=="return.level") se.theta <- sqrt(diag(var.theta))
	    out <- cbind(p - z.alpha * se.theta, p, p + z.alpha * se.theta)
	    rownames(out) <- par.name
	    conf.level <- paste(round((1 - alpha)*100, digits=2), "%", sep="")
	    colnames(out) <- c(paste(conf.level, " lower CI", sep=""), "Estimate", paste(conf.level, " upper CI", sep=""))
	    attr(out, "data.name") <- x$call
    	    attr(out, "method") <- method.name
    	    attr(out, "conf.level") <- (1 - alpha) * 100
    	    class(out) <- "ci"

    	    return(out)

	} else out <- c(p - z.alpha * sqrt(var.theta[ which.par ]), p, p + z.alpha * sqrt(var.theta[ which.par ])) # CJP2 - using var.theta

    } else if(method=="boot") {

	method.name <- "Parametric Bootstrap"

	if(verbose) cat("\n", "Simulating data from fitted model.  Size = ", R, "\n")
	if(const) {

	    if(is.null(x$blocks)) { # CJP

              Z <- rextRemes(x, n=R * x$n)
              Z <- matrix(Z, x$n, R)

            } else {

              Z <- rextRemes(x, n = round(R*x$npy*x$blocks$nBlocks))
              Z <- matrix(Z, round(x$npy*x$blocks$nBlocks), R)

            }

	} else Z <- rextRemes(x, n=R)
	if(verbose) cat("\n", "Simulated data found.\n")

	y <- datagrabber(x, response=FALSE)
	if(is.element(x$type, c("PP", "GP", "Exponential", "Beta", "Pareto"))) {

	    x2 <- datagrabber(x, cov.data=FALSE)
	    eid <- x2 > x$threshold
	    Z2 <- matrix(x$threshold, x$n, R)
	    Z2[eid,] <- Z[eid,]
	    Z <- Z2
	    lam <- mean(eid)

	} else {

	    eid <- !logical(x$n)
	    lam <- 1

	}
	ipar <- list()
	if(any(is.element(c("location","mu0"), theta.names))) {
	    if(is.element("location", theta.names)) ipar$location <- theta.hat["location"]
	    else {
		id <- substring(theta.names,1,2) == "mu"
		ipar$location <- theta.hat[id]
	    }
	}

	if(is.element("scale",theta.names)) ipar$scale <- theta.hat["scale"]
	else {
	    if(!x$par.models$log.scale) id <- substring(theta.names,1,3) == "sig"
	    else id <- substring(theta.names, 1, 3) == "phi"
	    ipar$scale <- theta.hat[id]
	}

	if(!is.element(x$type, c("Gumbel","Exponential"))) {
	    if(is.element("shape",theta.names)) ipar$shape <- theta.hat["shape"]
	    else {
		id <- substring(theta.names, 1, 2) == "xi"
		ipar$shape <- theta.hat[id]
	    }
	}

	bfun <- function(z, x, y, p, ipar, eid, rate) {
	    pm <- x$par.models
	    if(is.null(y)) fit <- fevd(x=z, threshold=x$threshold, location.fun=pm$location, scale.fun=pm$scale, shape.fun=pm$shape,
					    use.phi=pm$log.scale, type=x$type, method=x$method, initial=ipar, span=x$span,
					    time.units=x$time.units, period.basis=x$period.basis, optim.args=x$optim.args,
					    priorFun=x$priorFun, priorParams=x$priorParams, verbose=FALSE)
	    else fit <- fevd(x=z, data=y, threshold=x$threshold, location.fun=pm$location, scale.fun=pm$scale, shape.fun=pm$shape,
                                            use.phi=pm$log.scale, type=x$type, method=x$method, initial=ipar, span=x$span,
                                            time.units=x$time.units, period.basis=x$period.basis, optim.args=x$optim.args,
                                            priorFun=x$priorFun, priorParams=x$priorParams, verbose=FALSE)

	    fit$cov.data <- y
	    res <- distill(fit, cov=FALSE)

	    return(res)

	} # end of internal 'bfun' function.
	if(verbose) cat("\n", "Fitting model to simulated data sets (this may take a while!).")

	if(type=="parameter") {

	    sam <- apply(Z, 2, bfun, x=x, y=y, ipar=ipar)

	    if(tscale) {

                if(!const && !is.element("scale",theta.names) && !is.element("shape",theta.names)) stop("ci: invalid argument configurations.")
                
                if(!is.element(x$type, c("GP","Beta","Pareto"))) stop("ci: invalid argument configurations.")
                sam["scale",] <- sam["scale",] - sam["shape",] * x$threshold
		theta.hat["scale"] <- theta.hat["scale"] - theta.hat["shape"] * x$threshold
                theta.names[theta.names == "scale"] <- "tscale"
		rownames(sam) <- theta.names
		names(theta.hat) <- theta.names

	    } # end of if 'tscale' stmts.

	    sam <- sam[which.par,]

	    if(return.samples) return(t(sam))

        } else if(type=="return.level") {

	    pars <- apply(Z, 2, bfun, x=x, y=y, ipar=ipar)[1:np,]
	    th.est <- numeric(3)

	    if(is.element(x$type, c("PP","GEV","Gumbel","Weibull","Frechet"))) {
		loc <- pars["location",]
		th.est[1] <- theta.hat["location"]
	    } else loc <- rep(0, R)

	    scale <- pars["scale",]
	    th.est[2] <- theta.hat["scale"]

	    if(!is.element(x$type, c("Gumbel","Exponential"))) {
		shape <- pars["shape",]
		th.est[3] <- theta.hat["shape"]
	    } else {
		shape <- rep(1e-10, R)
		th.est[3] <- 1e-8
	    }

	    if(return.samples) out <- t(pars)

	    th <- rbind(loc, scale, shape)
	    rlfun <- function(theta, p, u, type, npy, rate)  rlevd(period=p, loc=theta[1], scale=theta[2], shape=theta[3], threshold=u, type=type, npy=npy, rate=rate)
	    if(x$type=="PP") mod <- "GEV"
	    else mod <- x$type
	    sam <- apply(th, 2, rlfun, p=return.period, u=x$threshold, type=mod, npy=x$npy, rate=lam)
	    if(is.matrix(sam)) rownames(sam) <- paste(rownames(sam), "-", x$period.basis, sep="")
	    else sammy.name <- paste(return.period, "-", x$period.basis, sep="")

	    if(return.samples) {

		if(is.matrix(sam)) out <- cbind(pars, t(sam))
		else {
		    onames <- colnames(out)
		    out <- cbind(out, sam)
		    colnames(out) <- c(onames, sammy.name)
		}
		return(out)
	    }

	    theta.hat <- rlevd(period=return.period, loc=th.est[1], scale=th.est[2], shape=th.est[3], threshold=x$threshold, type=x$type, npy=x$npy, rate=lam)

        } else stop("ci: invalid type argument.  Must be return.level or parameter.")

	if(is.matrix(sam)) {

	    out <- apply(sam, 1, quantile, probs=c(alpha/2, 1 - alpha/2))
            out.names <- rownames(out)
            out <- rbind(out[1,], theta.hat, out[2,])
            rownames(out) <- c(out.names[1], "Estimate", out.names[2])
	    colnames(out) <- rownames(sam)
            out <- t(out)
            attr(out, "data.name") <- x$call
            attr(out, "method") <- method.name
            attr(out, "conf.level") <- (1 - alpha) * 100
	    attr(out, "R") <- R
            class(out) <- "ci"

            return(out)

	} else {

	    out <- quantile(sam, probs=c(alpha/2, 1 - alpha/2))
            out <- c(out[1], mean(sam), out[2])
	    attr(out, "R") <- R

	}
	if(verbose) cat("\n", "Finished fitting model to simulated data.\n")

    } else if( method == "proflik" ) {

      if(x$type == "PP" && !is.null(x$blocks)) stop("ci: cannot do profile likelihood with blocks.") # CJP
      
	if(tscale) stop("ci: invalid argument configurations.")
	
	if( type == "parameter" && length( which.par ) > 1 ) stop("ci: can only do one parameter at a time with profile likelihood method.")
	else if( type == "return.level" && length( return.period ) > 1 ) stop("ci: can only do one return level at a time with profile likelihood method.")
	method.name <- "Profile Likelihood"
	if(verbose) {
	    if(x$type != "PP") cat("\n", "Calculating profile likelihood.  This may take a few moments.\n")
	    else cat("\n", "Calculating profile likelihood.  This may take several moments.\n")
	}

	if(is.null(xrange)) {

	    hold2 <- c(ci(x, alpha=alpha, method="normal", type=type, return.period=return.period, which.par=which.par))[c(1,3)]
	    # if(!any(is.na(hold2))) xrange <- c(p - hold2[1], p + hold2[2])
	    if(!any(is.na(hold2))) xrange <- range(c(hold2, log2(hold2), 4 * hold2, hold2 - 4 * hold2, hold2 + 4 * hold2), finite = TRUE)
	    else if(!is.na(hold2[2])) xrange <- range(c(p - 2 * abs(log2(abs(p))), hold2[2], 4 * hold2[2], -4 * hold2[2], log2(p)), finite = TRUE)
	    else if(!is.na(hold2[1])) xrange <- range(c(p - 2 * abs(log2(abs(p))), hold2[1], 4 * hold2[1], -4 * hold2[1], log2(p)), finite = TRUE)
	    else if(all(is.na(hold2))) xrange <- c(p - 2 * abs(log2(abs(p))), p + 2 * abs(log2(abs(p))))
	    # else if(!is.na(hold2[2])) xrange <- c(p - 2 * abs(log2(abs(p))), p + hold2[2])
	    # else if(!is.na(hold2[1])) xrange <- c(p - hold2[1], p + 2 * abs(log2(abs(p))))
	    # else if(all(is.na(hold2))) xrange <- c(p - 2 * abs(log2(abs(p))), p + 2 * abs(log2(abs(p))))
	    if(verbose) cat("\n", "Using a range of ", xrange[1], " to ", xrange[2], "\n")

	}

        if(is.null(x$blocks)) {

	    if(!is.null(xrange)) hold <- profliker(x, type=type, xrange=xrange, return.period=return.period, which.par=which.par, nint=nint, plot=verbose, ...)
	    else hold <- profliker(x, type=type, return.period=return.period, which.par=which.par, nint=nint, plot=verbose, ...)

        } else stop("Sorry: profile likelihood with blocks is not supported.")
	# CJP: I haven't figured out if we can implement the profile likelihood with the blocks approach

	ma <- -x$results$value
	crit <- ma - 0.5 * qchisq(1 - alpha, 1)

	if(verbose) {

	    cat("\n", "Profile likelihood has been calculated.  Now, trying to find where it crosses the critical value = ", crit, "\n")
	    abline(h=crit, col="blue")

	}

	crit2 <- ma - 0.5 * qchisq((1 - alpha)+abs(log2(1-alpha))/2, 1)

	id <- hold > crit2

	z <- seq(xrange[1], xrange[2], length=length(hold))
	z <- z[id]
	parlik <- hold[id]
	smth <- spline(z, parlik, n=200)
	ind <- smth$y > crit 
	out <- range(smth$x[ind])
	if(verbose) abline(v=out, lty=2, col="darkblue", lwd=2)
	out <- c(out[1], p, out[2])

    } else stop("ci: invalid method argument.")

    conf.level <- paste(round((1 - alpha)*100, digits=2), "%", sep="")
    names(out) <- c(paste(conf.level, " lower CI", sep=""), par.name, paste(conf.level, " upper CI", sep=""))
    attr(out, "data.name") <- x$call
    attr(out, "method") <- method.name
    attr(out, "conf.level") <- (1 - alpha) * 100
    class(out) <- "ci"

    return(out)

} # end of 'ci.fevd.mle' function.

return.level.ns.fevd.bayesian <- function(x, return.period = 100, ..., burn.in = 499, FUN = "mean", do.ci = FALSE, verbose = FALSE, qcov = NULL, qcov.base = NULL) {

    if(do.ci) res <- ci.rl.ns.fevd.bayesian(x = x, return.period = return.period, burn.in = burn.in, FUN = FUN, qcov = qcov, qcov.base = qcov.base, verbose = verbose, ...)
    else res <- erlevd(x, period = return.period)

    attr(res, "return.period") <- return.period
    attr(res, "data.name") <- x$data.name 
    attr(res, "fit.call") <- x$call
    attr(res, "call") <- match.call()
    attr(res, "fit.type") <- x$type
    attr(res, "data.assumption") <- "non-stationary"
    attr(res, "period") <- x$period.basis
    attr(res, "units") <- x$units
    attr(res, "class") <- "return.level"

    if(!do.ci) {

	attr(res, "conf.level") <- NULL
	class(res) <- "return.level"

    } else class(res) <- c("return.level", "ci")

    return(res)

} # end of 'return.level.ns.fevd.bayesian'

ci.rl.ns.fevd.bayesian <- function(x, alpha = 0.05, return.period = 100, FUN = "mean", burn.in = 499, ..., qcov = NULL, qcov.base = NULL, verbose = FALSE) {


    if(verbose) begin.tiid <- Sys.time()

    if(length(return.period) > 1) stop("ci.rl.ns.fevd.bayesian: return.period must have length 1.")

    conf.level <- paste(round((1 - alpha) * 100, digits = 2), "%", sep = "")

    np <- dim(x$results)[2] - 1
    if(burn.in > 0) p <- x$results[-(1:burn.in), 1:np]
    else p <- x$results[, 1:np]
    pnames <- colnames(p)

    ni <- dim(p)[1]

    if(is.null(qcov) && !is.null(qcov.base)) {

	warning("ci.rl.ns.fevd.bayesian: qcov must be supplied if qcov.base is supplied.  Continuing as if qcov.base were qcov, and qcov.base were NULL.")
	qcov <- qcov.base
	qcov.base <- NULL

    }

    rlfun <- function(th, pd, type, npy, rate) {

            return(rlevd(period = pd, loc = th[1], scale = th[2], shape = th[3], threshold = th[4], type = type, npy = npy, rate = rate))

    } # end of internal 'rlfun' function.

    if(!is.null(qcov)) {

	if(verbose) cat("\n", "Calculating ", return.period, "-year effective return levels based on qcov.\n")

	nq <- dim(qcov)[1]


	## Setting up parameter matrices so that they have the same number of rows as the MCMC (less burn.in)
	## and columns equal to the number of rows in the 'qcov' matrix.
	##

	# Set up location matrix
	if(is.element("location", pnames)) {

	    loc <- matrix(p[, "location"], ni, nq)
	    nloc <- 1

        } else if(is.element("mu", substring(pnames, 1, 2))) {

	    nloc <- sum(substring(pnames, 1, 2) == "mu")
	    loc <- matrix(NA, ni, nq)

	    for(i in 1:nq) loc[,i] <- rowSums(p[, 1:nloc, drop = FALSE] * matrix(qcov[i, 1:nloc], ni, nloc, byrow = TRUE))


        } else {

	    loc <- matrix(0, ni, nq)
	    nloc <- 0

        } # end of setting up location matrix stmts.


	# Set up scale matrix.
	if(is.element("log.scale", pnames)) {

	    nsc <- 1
	    scale <- matrix(exp(p[, "log.scale"]), ni, nq)

	} else if(is.element("scale", pnames)) {

	    nsc <- 1
	    scale <- matrix(p[, "scale"], ni, nq)

	} else if(is.element("sig", substring(pnames, 1, 3))) {

	    nsc <- sum(substring(pnames, 1, 3) == "sig")
	    scale <- matrix(NA, ni, nq)

	    for(i in 1:nq) scale[,i] <- rowSums(p[, (nloc + 1):(nloc + nsc), drop = FALSE] * matrix(qcov[i, (nloc + 1):(nloc + nsc)], ni, nsc, byrow = TRUE))

	} else if(is.element("phi", substring(pnames, 1, 3))) {

	    nsc <- sum(substring(pnames, 1, 3) == "phi")
	    scale <- matrix(NA, ni, nq)

	    for(i in 1:nq) scale[,i] <- rowSums(p[, (nloc + 1):(nloc + nsc), drop = FALSE] * matrix(qcov[i, (nloc + 1):(nloc + nsc)], ni, nsc, byrow = TRUE))

	    scale <- exp(scale)

	}

	# Set up shape matrix
	if(is.element("shape", pnames)) {

	    nsh <- 1
	    shape <- matrix(p[, "shape"], ni, nq)

	} else if(is.element("xi", substring(pnames, 1, 2))) {

	    nsh <- sum(substring(pnames, 1, 2) == "xi")
	    shape <- matrix(NA, ni, nq)

	    for(i in 1:nq) rowSums(p[, (nloc + nsc + 1):(nloc + nsc + nsh), drop = FALSE] * matrix(qcov[i, (nloc + nsc + 1):(nloc + nsc + nsh)], ni, nsh, byrow = TRUE))

	} else {

	    nsh <- 0
	    shape <- matrix(0, ni, nq)

	}

	u <- matrix(qcov[, "threshold"], ni, nq, byrow = TRUE)

	rlmat <- matrix(NA, ni, nq)

	for(i in 1:nq) rlmat[,i] <- apply(cbind(loc[,i], scale[,i], shape[,i], u[,i]), 1, rlfun, pd = return.period, type = x$type, npy = x$npy, rate = x$rate)

	if(!is.null(qcov.base)) {

	    if(dim(qcov.base)[1] != nq) stop("ci.rl.ns.fevd.bayesian: qcov.base must have same number of rows as qcov.")

	    # Set up location matrix
            if(is.element("location", pnames)) {

            loc <- matrix(p[, "location"], ni, nq)
            nloc <- 1

            } else if(is.element("mu", substring(pnames, 1, 2))) {

                nloc <- sum(substring(pnames, 1, 2) == "mu")
                loc <- matrix(NA, ni, nq)

                for(i in 1:nq) loc[,i] <- rowSums(p[, 1:nloc, drop = FALSE] * matrix(qcov.base[i, 1:nloc], ni, nloc, byrow = TRUE))


            } else {

                loc <- matrix(0, ni, nq)
                nloc <- 0

            } # end of setting up location matrix stmts.


            # Set up scale matrix.
            if(is.element("log.scale", pnames)) {

                nsc <- 1
                scale <- matrix(exp(p[, "log.scale"]), ni, nq)

            } else if(is.element("scale", pnames)) {

                nsc <- 1
                scale <- matrix(p[, "scale"], ni, nq)

            } else if(is.element("sig", substring(pnames, 1, 3))) {

                nsc <- sum(substring(pnames, 1, 3) == "sig")
                scale <- matrix(NA, ni, nq)

                for(i in 1:nq) scale[,i] <- rowSums(p[, (nloc + 1):(nloc + nsc), drop = FALSE] * matrix(qcov.base[i, (nloc + 1):(nloc + nsc)], ni, nsc, byrow = TRUE))

            } else if(is.element("phi", substring(pnames, 1, 3))) {

                nsc <- sum(substring(pnames, 1, 3) == "phi")
                scale <- matrix(NA, ni, nq)

                for(i in 1:nq) scale[,i] <- rowSums(p[, (nloc + 1):(nloc + nsc), drop = FALSE] * matrix(qcov.base[i, (nloc + 1):(nloc + nsc)], ni, nsc, byrow = TRUE))

                scale <- exp(scale)

            }

            # Set up shape matrix
            if(is.element("shape", pnames)) {

                nsh <- 1
                shape <- matrix(p[, "shape"], ni, nq)

            } else if(is.element("xi", substring(pnames, 1, 2))) {

                nsh <- sum(substring(pnames, 1, 2) == "xi")
                shape <- matrix(NA, ni, nq)

                for(i in 1:nq) rowSums(p[, (nloc + nsc + 1):(nloc + nsc + nsh), drop = FALSE] * matrix(qcov.base[i, (nloc + nsc + 1):(nloc + nsc + nsh)], ni, nsh, byrow = TRUE))

            } else {

                nsh <- 0
                shape <- matrix(0, ni, nq)

            }

            u <- matrix(qcov[, "threshold"], ni, nq, byrow = TRUE)

            rlmat2 <- matrix(NA, ni, nq)

	    for(i in 1:nq) rlmat2[,i] <- apply(cbind(loc[,i], scale[,i], shape[,i], u[,i]), 1, rlfun, pd = return.period, type = x$type, npy = x$npy, rate = x$rate)

	    res <- rlmat - rlmat2
	

	} else {

	    res <- rlmat

	} # end of if 'qcov.base' stmts.

    } else {

	if(verbose) cat("\n", "Calculating ", return.period, "-year effective return levels based on data covariates.\n")

	# Cannot use 'findpars' because we want parameters for every value of covariate at every iteration of MCMC chain.
	designs <- setup.design(x)

	X.loc <- designs$X.loc
	if(!is.null(X.loc)) nloc <- ncol(X.loc)
	else nloc <- 0

	X.sc <- designs$X.sc
	nsc <- ncol(X.sc)

	X.sh <- designs$X.sh
	if(!is.null(X.sh)) nsh <- ncol(X.sh)
	else nsh <- 0

	# loc <- scale <- shape <- res <- matrix(NA, x$n, ni)

	if(is.null(X.loc)) loc <- matrix(0, x$n, ni)
	if(is.null(X.sh)) shape <- matrix(0, x$n, ni)

	if(is.null(x$threshold)) u <- matrix(0, x$n, ni)
	else u <- matrix(x$threshold, x$n, ni)


	parfinder <- function(z, y) {

            return(rowSums(t(z * t(y)), na.rm = TRUE))

        } # end of internal 'parfinder' function.


	if(verbose) cat("\n", "Finding effective parameters for each iteration of MCMC sample.\n")

	if(is.null(X.loc)) loc <- matrix(0, ni, x$n)
        else loc <- apply(p[, 1:nloc, drop = FALSE], 1, parfinder, y = X.loc)

        scale <- apply(p[, (nloc+1):(nloc+nsc), drop = FALSE], 1, parfinder, y = X.sc)
        if(is.element("log.scale", pnames) || is.element("phi", substring(pnames, 1, 3))) scale <- exp(scale)

        if(is.null(X.sh)) shape <- matrix(0, ni, x$n)
        else shape <- apply(p[, (nloc+nsc+1):np, drop = FALSE], 1, parfinder, y = X.sh)

	do.rl <- function(id, z, u, pd, type, npy, rate, verbose) {

            if(verbose && id <= 5) cat(id, " ")
	    else if(verbose && (id %% 100 == 0)) cat(id, " ")

            theta <- cbind(z$loc[,id], z$scale[,id], z$shape[,id], z$threshold[,id])

	    id <- theta[,3] == 0
	    res <- numeric(dim(theta)[1]) + NA

            if(is.element(type, c("GEV","PP", "Gumbel"))) {

                p <- 1 - 1/pd
                if(any(id)) res[id] <- theta[id,1] - theta[id,2] * log(-log(p))
                if(any(!id)) res[!id] <- theta[!id,1] + theta[!id,2] * ((-log(p))^(-theta[!id,3]) - 1)/theta[!id,3]

            } else if(is.element(type, c("GP", "Exponential"))) {

		m <- pd * npy * rate
		if(any(id)) res[id] <- theta[id,4] + theta[id,2] * log(m)
		if(any(!id)) res[!id] <- theta[!id,4] + (theta[!id,2]/theta[!id,3]) * (m^(theta[!id,3]) - 1)

	    }

            return(res)

        } # end of internal 'do.rl' function.

	hold <- list(loc = loc, scale = scale, shape = shape, threshold = u)
	if(verbose) cat("\n", "Calculating the return levels for ", ni, " samples after burn in period.\n")
        res <- t(apply(matrix(1:ni, ncol = 1), 1, do.rl, z = hold, pd = 100, type = x$type, npy = x$npy, rate = x$rate, verbose = verbose))

    } # end of if else 'qcov' stmts.

    out <- cbind(c(apply(res, 2, quantile, probs = alpha / 2)), c(apply(res, 2, FUN, ...)), c(apply(res, 2, quantile, probs = 1 - alpha / 2)))
    colnames(out) <- c(paste(conf.level, " lower CI", sep = ""), "Estimate", paste(conf.level, " upper CI", sep = ""))

    attr(out, "conf.level") <- (1 - alpha) * 100
    class(out) <- "ci"

    attr(out, "return.period") <- return.period
    attr(out, "data.name") <- x$data.name
    attr(out, "fit.call") <- x$call
    attr(out, "call") <- match.call()
    attr(out, "fit.type") <- x$type
    attr(out, "data.assumption") <- "non-stationary"
    attr(out, "period") <- x$period.basis
    attr(out, "units") <- x$units
    attr(out, "class") <- "ci"

    if(verbose) {

	cat("\n", "Return levels and CIs estimated.\n")
	print(Sys.time() - begin.tiid)

    }

    return(out)

} # end of 'ci.rl.ns.fevd.bayesian' function.

ci.rl.ns.fevd.mle <- function(x, alpha = 0.05, return.period = 100, method = c("normal"), verbose = FALSE, qcov = NULL, qcov.base = NULL, ...) {

    method <- tolower(method)
    method <- match.arg(method)

    # if(method != "normal") stop("ci.rl.ns.fevd.mle: Currently only normal approximation method available for non-stationary MLE models.")

    par.name <- paste(return.period, "-", x$period.basis, " return level", sep = "")
    if (verbose) {

            cat("\n", "Preparing to calculate ", (1 - alpha) *
                100, "% CI for ", paste(return.period, "-", x$period.basis,
                " return level", sep = ""), "\n")

            cat("\n", "Model is non-stationary.\n")

            # if(method == "normal")
	    cat("\n", "Using Normal Approximation Method.\n")
            # else if(method == "boot") cat("\n", "Using Parametric Boot strap method with ", R, " replicate samples.\n")

    }

    if(method == "normal") method.name <- "Normal Approx."
    # else if(method == "boot") method.name <- "Parametric Bootstrap"

    if(method == "normal") {

        res <- return.level.ns.fevd.mle(x = x, return.period = return.period, ..., do.ci = FALSE, verbose = verbose, qcov = qcov, qcov.base = qcov.base)

        z.alpha <- qnorm(alpha/2, lower.tail = FALSE)
        cov.theta <- parcov.fevd(x)

        if (is.null(cov.theta))
            stop("ci: Sorry, unable to calculate the parameter covariance matrix.  Maybe try a different method.")
        var.theta <- diag(cov.theta)
        if (any(var.theta < 0))
            stop("ci: negative Std. Err. estimates obtained.  Not trusting any of them.")
        grads <- t(rlgrad.fevd(x, period = return.period, qcov = qcov,
            qcov.base = qcov.base))
        se.theta <- sqrt(diag(t(grads) %*% cov.theta %*% grads))
        out <- cbind(c(res) - z.alpha * se.theta, c(res), c(res) + z.alpha * se.theta, se.theta)
        if (length(return.period) > 1)
            rownames(out) <- par.name
        else rownames(out) <- NULL
        conf.level <- paste(round((1 - alpha) * 100, digits = 2),
            "%", sep = "")
        colnames(out) <- c(paste(conf.level, " lower CI", sep = ""),
            "Estimate", paste(conf.level, " upper CI", sep = ""),
            "Standard Error")
    } else if(method == "boot") {

	stop("ci.rl.ns.fevd.mle: Sorry, this functionality has not yet been added.")

    }

    attr(out, "data.name") <- x$call
    attr(out, "method") <- method.name
    attr(out, "conf.level") <- (1 - alpha) * 100
    class(out) <- "ci"

    return(out)

} # end of 'ci.rl.ns.fevd.mle' function.


rlgrad.fevd <- function(x, period=100, qcov=NULL, qcov.base=NULL) {  # CJP2; many changes here for nonstationary models
  # qcov.base is for when one wants to work with the difference in return levels for different covariate sets in a nonstationary model

  type <- tolower(x$type)
  if(!is.element(x$method, c("MLE","GMLE"))) stop("rlgrad.fevd: Estimation method must be MLE/GMLE.")
  p <- x$results$par
  if(is.element("log.scale",names(p))) {
    id <- names(p) == "log.scale"
    p[id] <- exp(p[id])
    names(p)[id] <- "scale"
  } 

  if(is.element("shape",names(p))) {
    if(p["shape"] == 0) {
      if(is.element(type, c("gev","pp","gumbel"))) type <- "gumbel"
      else if(is.element(type, c("gp","exponential"))) type <- "exponential"
      else stop("rlgrad.fevd: invalid type for the shape parameter.")
    }
  }
  
  if(!is.fixedfevd(x)) { # CJP2: this whole block of code
    if(is.null(qcov)) stop("rlgrad.fevd: qcov required for nonstationary models.")
    if(!is.matrix(qcov)) qcov <- matrix(qcov, nrow = 1)
    if(!is.null(qcov.base)) {
      if(!is.matrix(qcov.base)) qcov.base <- matrix(qcov.base, nrow = 1)
      if(nrow(qcov) != nrow(qcov.base) ||
         ncol(qcov) != ncol(qcov.base))
        stop("rlgrad.fevd: qcov and qcov.base must have the same number of covariates and values.")
    }
    
    if(length(period) > 1 && nrow(qcov) > 1)
      stop("rlgrad.fevd: Cannot compute gradient for multiple return periods and multiple covariate values simultaneously.")

    loc.id <- 1:x$results$num.pars$location
    sc.id <- (1 + x$results$num.pars$location) : (x$results$num.pars$location + x$results$num.pars$scale)
    sh.id <- (1 + x$results$num.pars$location + x$results$num.pars$scale) : (x$results$num.pars$location + x$results$num.pars$scale + x$results$num.pars$shape)

    # loc <- qcov[ , loc.id, drop=FALSE] %*% p[loc.id]  # not needed
    sc <- c(qcov[ , sc.id, drop=FALSE] %*% p[sc.id])
    sh <- c(qcov[ , sh.id, drop=FALSE] %*% p[sh.id])

    if(!is.null(qcov.base)) {
      # loc.base <- qcov.base[ , loc.id, drop=FALSE] %*% p[loc.id]  # not needed
      sc.base <- c(qcov.base[ , sc.id, drop=FALSE] %*% p[sc.id])
      sh.base <- c(qcov.base[ , sh.id, drop=FALSE] %*% p[sh.id])
    }

    # CJP2 - this next block of code deals with the only case that params are on log scale when considering return levels, namely when use.phi=TRUE for nonstationarity in the scale
    pnames <- names(p)
    if(is.element("phi", substring(pnames, 1, 3))) {
      sc <- exp(sc)
      phi.gradTerm <- sc
      if(!is.null(qcov.base)) {
        sc.base <- exp(sc.base)
        phi.gradTerm.base <- sc.base
      }
    } else{
      phi.gradTerm <- 1
      if(!is.null(qcov.base)) phi.gradTerm.base <- 1
    }
    
    if(is.element(type, c("pp","gev","weibull","frechet"))) {
      yp <- -log(1 - 1/period)
      if(length(period) == 1) {
        res <- cbind(qcov[ , loc.id, drop = FALSE], -phi.gradTerm/sh * (1-yp^(-sh)) * qcov[ , sc.id, drop = FALSE], 
                     (sc * sh^(-2) * (1 - yp^(-sh)) - sc/sh * yp^(-sh) * log(yp)) * qcov[ , sh.id, drop = FALSE])
        if(!is.null(qcov.base))
          res <- res - cbind(qcov.base[ , loc.id, drop = FALSE], -phi.gradTerm.base/sh.base * (1-yp^(-sh.base)) * qcov.base[ , sc.id, drop = FALSE], (sc.base * sh.base^(-2) * (1 - yp^(-sh.base)) - sc.base/sh.base * yp^(-sh.base) * log(yp)) * qcov.base[ , sh.id, drop = FALSE])
      } else {
        res <- cbind(outer(rep(1, length(period)), qcov[1, loc.id]), outer(-phi.gradTerm/sh * (1-yp^(-sh)), qcov[1, sc.id]), outer(sc * sh^(-2) * (1 - yp^(-sh)) - sc/sh * yp^(-sh) * log(yp), qcov[1, sh.id]))
        if(!is.null(qcov.base))
          res <- res - cbind(outer(rep(1, length(period)), qcov.base[1, loc.id]), outer(-phi.gradTerm.base/sh.base * (1-yp^(-sh.base)), qcov.base[1, sc.id]), outer(sc.base * sh.base^(-2) * (1 - yp^(-sh.base)) - sc.base/sh.base * yp^(-sh.base) * log(yp), qcov.base[1, sh.id]))
      }
    } else if(type=="gumbel") {
      yp <- -log(1 - 1/period)
      if(length(period) == 1) {
        res <- cbind(qcov[ , loc.id, drop = FALSE], -log(yp)*phi.gradTerm * qcov[ , sc.id, drop = FALSE])
        if(!is.null(qcov.base))
          res <- res - cbind(qcov.base[ , loc.id, drop = FALSE], -log(yp)*phi.gradTerm.base * qcov.base[ , sc.id, drop = FALSE])
      } else {
        res <- cbind(outer(rep(1, length(period)), qcov[1, loc.id]), outer(-log(yp)*phi.gradTerm, qcov[1, loc.id]))
        if(!is.null(qcov.base))
          res <- res - cbind(outer(rep(1, length(period)), qcov.base[1, loc.id]), outer(-log(yp)*phi.gradTerm.base, qcov.base[1, loc.id]))          
      }
    } else stop("rlgrad.fevd: not implemented for nonstationary models for GP, beta, pareto, exponential models.")
  } else { # fixed.fevd  
    if(is.element(type, c("pp","gev","weibull","frechet"))) {
      yp <- -log(1 - 1/period)
      res <- cbind(1,
                   (-1/p["shape"]) * (1 - yp^(-p["shape"])), 
                   p["scale"] * (p["shape"])^(-2) * (1 - yp^(-p["shape"])) - (p["scale"]/p["shape"]) * yp^(-p["shape"]) * log(yp))
    } else if(type=="gumbel") {
      yp <- -log(1 - 1/period)
      res <- cbind(1, -log(yp)) 
    } else if(is.element(type, c("gp","beta","pareto"))) {
      lam <- mean(c(datagrabber(x)[,1]) > x$threshold)
      m <- period * x$npy
      mlam <- m * lam
      res <- cbind(p["scale"] * m^(-p["shape"]) * lam^(-p["shape"] - 1),
                   (p["shape"])^(-1) * ((mlam)^(p["shape"]) - 1),  
                   -p["scale"] * (p["shape"])^(-2)*((mlam)^(p["shape"]) - 1) + (p["scale"]/p["shape"]) * (mlam)^(p["shape"]) * log(mlam))
    } else if(type=="exponential") {
      lam <- mean(c(datagrabber(x)[,1]) > x$threshold)
      m <- period * x$npy 
      mlam <- m * lam
      res <- cbind(p["scale"]/mlam, log(mlam)) 
    } else stop("rlgrad.fevd: Hmmm.  You should not be getting this error message.  Something is horribly wrong.")
  } # end of if else use actual gradients or finite differences stmts.
  return(res)
} # end of 'rlgrad.fevd' function.

return.level <- function(x, return.period=c(2, 20, 100), ...) {
    UseMethod("return.level", x)
} # end of 'return.level' function.

return.level.fevd <- function(x, return.period=c(2, 20, 100), ...) {
  newcl <- tolower(x$method)
  if(newcl == "gmle") newcl <- "mle"
  class(x) <- paste("fevd.", newcl, sep="")
  UseMethod("return.level", x)
} # end of 'return.level.fevd' function.

return.level.fevd.lmoments <- function(x, return.period=c(2, 20, 100), ..., do.ci = FALSE) {
    model <- x$type

    p <- x$results
    pnames <- names(p)

    if(model=="PP") mod2 <- "GEV"
    else mod2 <- model

    if(!do.ci) {
	if(all(is.element(c("location","shape"),pnames))) {
            res <- rlevd(return.period, loc=p["location"], scale=p["scale"], shape=p["shape"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)
        } else if(is.element("shape", pnames)) {
            res <- rlevd(return.period, scale=p["scale"], shape=p["shape"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)
        } else if(is.element("location", pnames)) {
            res <- rlevd(return.period, loc=p["location"], scale=p["scale"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)
        } else res <- rlevd(return.period, scale=p["scale"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)
        attr(res, "return.period") <- return.period
        attr(res, "data.name") <- x$data.name
        attr(res, "fit.call") <- x$call
        attr(res, "call") <- match.call()
        attr(res, "fit.type") <- x$type
        attr(res, "data.assumption") <- "stationary"
        attr(res, "period") <- x$period.basis
        attr(res, "units") <- x$units
        attr(res, "class") <- "return.level"
    } else if(do.ci) res <- ci(x, return.period=return.period, ...)

    return(res)
} # end of 'return.level.fevd.lmoments' function.


return.level.fevd.bayesian <- function(x, return.period = c(2, 20, 100), ..., do.ci = FALSE, burn.in = 499, FUN = "mean", qcov = NULL, qcov.base = NULL) {

    model <- x$type
    if(model=="PP") mod2 <- "GEV"
    else mod2 <- model

    # tform <- !is.fixedfevd(x)  # CJP2

    if(model=="PP") mod2 <- "GEV"
    else mod2 <- model

    f <- match.fun(FUN)

    p <- x$results
    np <- dim(p)[2] - 1
    p <- p[,1:np]
    pnames <- colnames(p)

    if(FUN=="mean") p <- colMeans(p)
    else p <- apply(p, 2, f)

    if(is.fixedfevd(x)) {  # CJP2

        if(is.element("log.scale",pnames)) {
            p["log.scale"] <- exp(p["log.scale"])
            pnames[pnames == "log.scale"] <- "scale"
            names(p) <- pnames
        }

        if(!do.ci) {

            if(all(is.element(c("location","shape"),pnames))) {

                res <- rlevd(return.period, loc=p["location"], scale=p["scale"], shape=p["shape"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)

            } else if(is.element("shape", pnames)) {

                res <- rlevd(return.period, scale=p["scale"], shape=p["shape"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)

            } else if(is.element("location", pnames)) {

                res <- rlevd(return.period, loc=p["location"], scale=p["scale"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)

            } else res <- rlevd(return.period, scale=p["scale"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)

                attr(res, "return.period") <- return.period
                attr(res, "data.name") <- x$data.name
                attr(res, "fit.call") <- x$call
                attr(res, "call") <- match.call()
                attr(res, "fit.type") <- x$type
                attr(res, "data.assumption") <- "stationary"
                attr(res, "period") <- x$period.basis
                attr(res, "units") <- x$units
                attr(res, "class") <- "return.level"

            } else if(do.ci) res <- ci(x, return.period=return.period, ...)

    } else {

	if(missing(return.period)) return.period <- 100
	res <- return.level.ns.fevd.bayesian(x = x, return.period = return.period, ..., do.ci = do.ci, qcov = qcov, qcov.base = qcov.base)
	if(do.ci) return(res)

        if(length(return.period)==1) res <- matrix(res, ncol=1)
        colnames(res) <- paste(return.period, "-", x$period.basis, " level", sep="")
        attr(res, "return.period") <- return.period
        attr(res, "data.name") <- x$data.name
        attr(res, "fit.call") <- x$call
        attr(res, "call") <- match.call()
        attr(res, "fit.type") <- x$type
        attr(res, "data.assumption") <- "non-stationary"
        attr(res, "period") <- x$period.basis
        attr(res, "units") <- x$units
        if(is.null(qcov)) attr(res, "qcov") <- x$data.name[2]
        else attr(res, "qcov") <- deparse(substitute(qcov))
        attr(res, "class") <- "return.level"
    }

    return(res)

} # end of 'return.level.fevd.bayesian' function.

return.level.ns.fevd.mle <-
function (x, return.period = c(2, 20, 100), ..., alpha = 0.05,
    method = c("normal"), do.ci = FALSE, verbose = FALSE, qcov = NULL,
    qcov.base = NULL)
{

    if(missing(return.period)) return.period <- 100

    if (do.ci && method != "normal")
        stop("return.level.ns.fevd.mle: only normal approximation CI calculations currently available for nonstationary return levels.")

    model <- x$type

    if(!(tolower(model) %in% c("pp", "gev", "weibull", "frechet", "gumbel")))
        stop("return.level.ns.fevd.mle: not implemented for GP, beta, pareto, exponential models.")

    if (do.ci && length(return.period) > 1 && nrow(qcov) > 1)
        stop("return.level.ns.fevd.mle:: Cannot calculate confidence intervals for multiple return periods and multiple covariate values simultaneously.")

    if(!do.ci) {

        if (model == "PP") {

            mod2 <- "GEV"

        } else mod2 <- model

        p <- x$results$par
        pnames <- names(p)

        if(is.fixedfevd(x)) stop("return.level.ns.fevd.mle: this function is for nonstationary models.")

        if(is.null(qcov)) stop("return.level.ns.fevd.mle: qcov required for this function.")

        if(!is.matrix(qcov)) qcov <- matrix(qcov, nrow = 1)

        if(!is.qcov(qcov)) qcov <- make.qcov(x = x, vals = qcov, nr = nrow(qcov))

        if(!is.null(qcov.base)) {

            if(!is.matrix(qcov.base)) qcov.base <- matrix(qcov.base, nrow = 1)

            if(nrow(qcov) != nrow(qcov.base) || ncol(qcov) != ncol(qcov.base))
                stop("return.level.ns.fevd.mle: qcov and qcov.base must have the same number of covariates and values.")

            if(!is.qcov(qcov.base)) qcov.base <- make.qcov(x = x, vals = qcov.base, nr = nrow(qcov.base))

        }

        loc.id <- 1:x$results$num.pars$location

        sc.id <- (1 + x$results$num.pars$location):(x$results$num.pars$location + x$results$num.pars$scale)

        sh.id <- (1 + x$results$num.pars$location + x$results$num.pars$scale):(x$results$num.pars$location +
            x$results$num.pars$scale + x$results$num.pars$shape)

        loc <- qcov[, loc.id, drop = FALSE] %*% p[loc.id]
        sc <- qcov[, sc.id, drop = FALSE] %*% p[sc.id]
        sh <- qcov[, sh.id, drop = FALSE] %*% p[sh.id]

        if(!is.null(qcov.base)) {

            loc.base <- qcov.base[, loc.id, drop = FALSE] %*% p[loc.id]
            sc.base <- qcov.base[, sc.id, drop = FALSE] %*% p[sc.id]
            sh.base <- qcov.base[, sh.id, drop = FALSE] %*% p[sh.id]

        }

        if(x$par.models$log.scale) {

            sc <- exp(sc)
            if(!is.null(qcov.base)) sc.base <- exp(sc.base)

        }

        theta <- cbind(qcov[, "threshold"], loc, sc, sh)

        rlfun2 <- function(th, pd, type, npy, rate) rlevd(pd, loc = th[2],
            scale = th[3], shape = th[4], threshold = th[1], type = type,
            npy = npy)

        res <- apply(theta, 1, rlfun2, pd = return.period, type = mod2, npy = x$npy)

        if (!is.null(qcov.base)) {

            theta.base <- cbind(qcov.base[, "threshold"], loc.base, sc.base, sh.base)

            res <- res - apply(theta.base, 1, rlfun2, pd = return.period, type = mod2, npy = x$npy)
        }

        res <- t(matrix(res, nrow = length(return.period)))
        colnames(res) <- paste(return.period, "-", x$period.basis,
            " level", sep = "")
        attr(res, "return.period") <- return.period
        attr(res, "data.name") <- x$data.name
        attr(res, "fit.call") <- x$call
        attr(res, "call") <- match.call()
        attr(res, "fit.type") <- x$type
        attr(res, "data.assumption") <- "non-stationary"
        attr(res, "period") <- x$period.basis
        attr(res, "units") <- x$units
        attr(res, "qcov") <- deparse(substitute(qcov))
        attr(res, "class") <- "return.level"
        if (!is.null(qcov.base)) {
            attr(res, "qcov.base") <- deparse(substitute(qcov.base))
            attr(res, "class") <- "return.level.diff"
        }
        return(res)

    } else {

        out <- ci.rl.ns.fevd.mle(x = x, alpha = alpha, return.period = return.period, qcov = qcov, qcov.base = qcov.base, ...)
        return(out)

    }

} # end of 'return.level.ns.fevd.mle' function.

return.level.fevd.mle <- function(x, return.period = c(2, 20, 100), ..., do.ci = FALSE, qcov = NULL, qcov.base = NULL) {

  model <- x$type

  # tform <- !is.fixedfevd(x) # removed by CJP2

  if(model=="PP") mod2 <- "GEV"
  else mod2 <- model
  
  p <- x$results$par
  pnames <- names(p)
  
  if(is.fixedfevd(x)) {  # CJP2
    
    if(is.element("log.scale",pnames)) {
      p["log.scale"] <- exp(p["log.scale"])
      pnames[pnames == "log.scale"] <- "scale"
      names(p) <- pnames
    }
    
    if(!do.ci) {
      if(all(is.element(c("location","shape"),pnames))) {
        res <- rlevd(return.period, loc=p["location"], scale=p["scale"], shape=p["shape"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)
      } else if(is.element("shape", pnames)) {
        res <- rlevd(return.period, scale=p["scale"], shape=p["shape"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)
      } else if(is.element("location", pnames)) {
        res <- rlevd(return.period, loc=p["location"], scale=p["scale"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)
      } else res <- rlevd(return.period, scale=p["scale"], threshold=x$threshold, type=mod2, npy=x$npy, rate=x$rate)
      attr(res, "return.period") <- return.period
      attr(res, "data.name") <- x$data.name
      attr(res, "fit.call") <- x$call
      attr(res, "call") <- match.call()
      attr(res, "fit.type") <- x$type
      attr(res, "data.assumption") <- "stationary"
      attr(res, "period") <- x$period.basis
      attr(res, "units") <- x$units
      attr(res, "class") <- "return.level"
    } else if(do.ci) res <- ci(x, return.period=return.period, ...)
    
  } else {

    if(missing(return.period)) return.period <- 100

    if(do.ci) {

	res <- ci.rl.ns.fevd.mle(x = x, return.period = return.period, ..., qcov = qcov, qcov.base = qcov.base)

    } else {

	if(is.null(qcov) && !is.null(qcov.base)) {
	    qcov <- qcov.base
	    qcov.base <- NULL
	    warning("return.level.fevd.mle: attempt to set qcov to null but not qcov.base.  Setting qcov to qcov.base and qcov.base to NULL.")
	}

        if(is.null(qcov)) {

          class(x) <- "fevd"
          rlfun <- function(p, x) return(erlevd(x=x, period=p))
          res <- apply(matrix(return.period, ncol=1), 1, rlfun, x=x)

        } else {
	    
	    res <- return(return.level.ns.fevd.mle(x = x, return.period = return.period, ..., do.ci = FALSE, qcov = qcov, qcov.base = qcov.base))
	    # return(return.level.ns.fevd.mle(x = x, return.period = return.period, ..., do.ci = do.ci, qcov = qcov, qcov.base = qcov.base))
# 
#           if(!is.matrix(qcov)) qcov <- matrix(qcov, nrow=1)
#           if(!is.qcov(qcov)) qcov <- make.qcov(x=x, vals=qcov, nr=nrow(qcov))
#       
#           nr <- nrow(qcov)
#       
# 	    if(is.element("location", pnames)) {
# 
# 	        nloc <- 1
#                 loc <- p["location"]
# 
# 	    } else if(is.element("mu", substring(pnames,1,2))) {
# 
#             id <- substring(pnames,1,2) == "mu"
#             nloc <- sum(id)
#             loc <- rowSums(matrix(p[1:nloc], nr, nloc, byrow=TRUE) * qcov[,1:nloc])
# 
#           } else loc <- nloc <- 0
#       
# 	if(is.element("scale", pnames)) {
# 
#             nsc <- 1
#             scale <- p["scale"]
# 
# 	} else if(is.element("phi", substring(pnames, 1, 3))) {
# 
#             id <- substring(pnames,1,3) == "phi"
#             nsc <- sum(id)
#             scale <- exp(rowSums(matrix(p[(nloc+1):(nloc+nsc)], nr, nsc, byrow=TRUE) * qcov[,(nloc+1):(nloc+nsc)]))
# 
# 	} else if(is.element("sig", substring(pnames, 1, 3))) {
# 
#         id <- substring(pnames,1,3) == "sig"
#         nsc <- sum(id)
#         scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], nr, nsc, byrow=TRUE) * qcov[,(nloc+1):(nloc+nsc)])
# 
# 	}
#       
# 	if(is.element("shape", pnames)) {
# 
#             nsh <- 1
#             shape <- p["shape"]
# 
# 	} else if(is.element("xi", substring(pnames,1,2))) {
# 
#             id <- substring(pnames,1,2) == "xi"
#             nsh <- sum(id)
#             shape <- rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], nr, nsh, byrow=TRUE) * qcov[,(nloc+nsc+1):(nloc+nsc+nsh)])
# 
# 	} else nsh <- shape <- 0
#       
# 	    theta <- cbind(qcov[,"threshold"], loc, scale, shape)
# 	    rlfun2 <- function(th, pd, type, npy, rate) rlevd(return.period, loc=th[2], scale=th[3], shape=th[4], threshold=th[1], type=type, npy=npy, rate=rate)
# 
# 	    res <- t(apply(theta, 1, rlfun2, pd=return.period, type=mod2, npy=x$npy, rate=x$rate))
# 
	} # end of if else 'qcov' argument is NULL stmts.
    
	if(length(return.period)==1) res <- matrix(res, ncol=1)
	colnames(res) <- paste(return.period, "-", x$period.basis, " level", sep="")
	attr(res, "return.period") <- return.period
	attr(res, "data.name") <- x$data.name
	attr(res, "fit.call") <- x$call
	attr(res, "call") <- match.call()
	attr(res, "fit.type") <- x$type
	attr(res, "data.assumption") <- "non-stationary"
	attr(res, "period") <- x$period.basis
	attr(res, "units") <- x$units
	if(is.null(qcov)) attr(res, "qcov") <- x$data.name[2]
	else attr(res, "qcov") <- deparse(substitute(qcov))
	attr(res, "class") <- "return.level"

	} # end of if else 'do.ci' stmts.

    } # end of if else fixed model stmts.
    
    return(res)

} # end of 'return.level.fevd.mle' function.

print.return.level <- function(x, ...) {

    tmp <- attributes(x)
    print(tmp$fit.call)
    print(tmp$call)

    if(!is.null(tmp$units)) print(paste(tmp$fit.type, " model fitted to ", tmp$data.name[1], " (", tmp$units, ")", sep=""))
    else cat("\n", tmp$fit.type, "model fitted to ", tmp$data.name, "\n")

    if(tmp$fit.type=="PP") print(paste("Return levels based on GEV equivalency (i.e., return levels are for block maxima, where the blocks are ", tmp$period, "s)", sep=""))
    cat("Data are assumed to be ", tmp$data.assumption, "\n")

    if(tmp$data.name[2] != "") {

	print(paste("Covariate data = ", tmp$data.name[2], sep=""))
	if(!is.null(tmp$qcov)) if(tmp$qcov != tmp$data.name[2]) print(paste("Covariate data used for effective return levels here = ", tmp$qcov, sep=""))

    }

    print(paste("Return Levels for period units in ", tmp$period, "s", sep=""))
    # cat("Return Period(s) = \n", tmp$return.period, "\n")
    if(is.null(tmp$dim)) {

        y <- c(x)
        names(y) <- paste(tmp$return.period, "-", tmp$period, " level", sep="")

    } else {

	y <- matrix(c(x), tmp$dim[1], tmp$dim[2])
	colnames(y) <- paste(tmp$return.period, "-", tmp$period, " level", sep="")

    }

    print(y)
    invisible()

} # end of 'print.return.level' function.

make.qcov <- function(x, vals, nr=1, ...) {

    if(x$method != "Bayesian") p <- x$results$par
    else p <- apply(x$results[, -dim(x$results)[2] ], 2, mean, na.rm=TRUE)


    np <- length(p)
    pnames <- names(p)

    if(is.element("mu", substring(pnames, 1, 2))) nloc <- sum(substring(pnames, 1, 2) == "mu")
    else if(is.element("location", pnames)) nloc <- 1
    else nloc <- 0

    if(is.element("phi", substring(pnames, 1, 3))) nsc <- sum(substring(pnames, 1, 3) == "phi")
    else if(is.element("sig", substring(pnames, 1, 3))) nsc <- sum(substring(pnames, 1, 3) == "sig")
    else nsc <- 1

    if(is.element("shape", pnames)) nsh <- 1
    else if(is.element("xi", substring(pnames, 1, 2))) nsh <- sum(substring(pnames, 1, 2) == "xi")
    else nsh <- 0

    if(missing(vals)) {
        out <- matrix(0, nr, np+1)
        out[,1] <- 1
        if(nloc > 0) out[,nloc+1] <- 1
        out[,nloc+nsc+1] <- 1
        if(!is.null(x$threshold)) {
	    if(length(x$threshold) >= nr) out[,np+1] <- x$threshold[1:nr]
	    else out[,np+1] <- x$threshold[1]
	}
        colnames(out) <- c(pnames, "threshold")
	return(out)
    }

    if(is.list(vals)) {

	vnames <- names(vals)
	if(is.null(vnames)) stop("make.qcov: vals must be a named list.")
	if(!all(is.element(vnames, c(pnames, "threshold")))) stop("make.qcov: names of vals list must match parameter names.")

	if(missing(nr)) {

	    nv <- lapply(vals, length)
	    if(length(unique(nv)) != 1) stop("make.qcov: Sorry, this function is limited.  Length of each component must be the same or nr must be specified.")

	    nr <- length(vals[[1]])
        }

	out <- matrix(NA, nrow=nr, ncol=np + 1)

	for(i in 1:np) {

	    if(is.element(pnames[i], vnames)) {
		id <- (1:length(vnames))[ vnames == pnames[i] ]
		out[,i] <- vals[[ id ]]
	    } else out[,i] <- 1

	} # end of for 'i' loop.

	if(is.element("threshold", vnames)) {

	    id <- (1:length(vnames))[ vnames == "threshold" ]
	    out[,np + 1] <- vals[[ id ]]

	} else if(!is.element(x$type, c("GEV", "Gumbel", "Frechet"))) out[,np + 1] <- x$threshold[1]
	

    } else if(is.matrix(vals)) out <- vals
    else if(is.numeric(vals)) out <- matrix(vals, nrow=nr, ...)
    else stop("make.qcov: invalid vals argument.")

    if(dim(out)[2]==np) {
	if(!is.null(x$threshold)) {
	    if(length(x$threshold) >= nr) out <- cbind(out, x$threshold[1:nr])
	    else out <- cbind(out, x$threshold[1])
	} else out <- cbind(out, 0)
	colnames(out) <- c(pnames, "threshold")
    } else if(dim(out)[2]==np+1) {
	if(is.null(colnames(out))) colnames(out) <- c(pnames, "threshold")
    } else stop("make.qcov: length of (or number of columns of) vals must be equal to number of model parameters, np, or np + 1.")

    # Some final checks.
    if(is.element("location", pnames)) {
	if(!all(out[,"location"]==1)) {
	    warning("make.qcov: invalid qcov values for fixed location parameter.  Re-setting to one.")
	    out[,"location"] <- 1
	}
    } else if(is.element("mu0", pnames)) {
	if(!all(out[,"mu0"] == 1)) {
	    warning("make.qcov: invalid qcov values for fixed mu0 parameter.  Re-setting to one.")
            out[,"mu0"] <- 1
	}
    }

    if(any(is.element(c("log.scale","scale"), pnames))) {
	id <- (pnames == "scale") | (pnames == "log.scale")
	id <- c(id, FALSE)
	if(!all(out[,id] == 1)) {
	    warning("make.qcov: invalid qcov values for fixed scale parameter.  Re-setting to one.")
            out[,id] <- 1
	}
    } else if(any(is.element(c("sig0","phi0"), pnames))) {
	id <- (pnames == "sig0") | (pnames == "phi0")
        id <- c(id, FALSE)
	if(!all(out[,id] == 1)) {
            warning("make.qcov: invalid qcov values for fixed scale intercept term.  Re-setting to one.")
            out[,id] <- 1
        }
    }

    if(is.element("shape", pnames)) {
	if(!all(out[,"shape"] == 1)) {
	    warning("make.qcov: invalid qcov values for fixed shape parameter.  Re-setting to one.")
            out[,"shape"] <- 1
	}
    } else if(is.element("xi0", pnames)) {
	if(!all(out[,"xi0"] == 1)) {
	    warning("make.qcov: invalid qcov values for fixed xi0 parameter.  Re-setting to one.")
            out[,"xi0"] <- 1
	}
    }

    return(out)
} # end of 'make.qcov' function.

is.qcov <- function(x) {
    if(!is.matrix(x)) return(FALSE)
    d <- dim(x)
    if(is.null(d)) return(FALSE)
    dnames <- colnames(x)
    if(is.null(dnames)) return(FALSE)
    if(dnames[d[2]] != "threshold") return(FALSE)
    if(d[2] > 1) return(TRUE)
    else return(FALSE)
} # end of 'is.qcov' function.

probprob.plot.evd <- function(xp, y, model, loc, scale, shape, u, tform=FALSE, eid, obj, ytrans=NULL, npy=NULL, ...) {
	args <- list(...)
        if(is.null(args$main))  m1 <- deparse(obj$call)

        if(!tform) {
            if(is.element(model, c("PP","GP","Beta","Pareto","Exponential"))) {
                yp <- pevd(sort(y[eid]), loc=loc, scale=scale, shape=shape, threshold=u, npy=npy, type=model)
            } else yp <- pevd(sort(y), loc=loc, scale=scale, shape=shape, threshold=u, npy=npy, type=model)
            if(is.null(args$main)) plot(xp, yp, main=m1, xlab="Empirical Probabilities", ylab="Model Probabilities", ...)
            else  plot(xp, yp, xlab="Empirical Probabilities", ylab="Model Probabilities", ...)
            abline(0,1)
        } else {
            # yp <-  pevd(sort(ytrans), loc=0, scale=1, shape=0, threshold=0, npy=npy, type=model)
            if(is.element(model, c("GEV","Gumbel","Weibull","Frechet"))) yp <- exp(-exp(-sort(ytrans)))
            else if(model=="PP") yp <- sort(ytrans)
            else yp <- 1 - exp(-sort(ytrans))
            if(is.null(args$main)) {
                plot(xp, yp, main=m1, xlab="Residual Empirical Probabilities", ylab="Residual Model Probabilities", ...)
            } else  plot(xp, yp, xlab="Residual Empirical Probabilities", ylab="Residual Model Probabilities", ...)
            abline(0,1)
        }
    invisible()
} # end of 'probprob.plot.evd' stmts.

quantquant.plot.evd <- function(x, xp, y, u, loc, scale, shape, tform=FALSE, eid, ytrans=NULL,
				    model= c("GEV", "GP", "PP", "Gumbel", "Frechet", "Weibull", "Exponential", "Beta", "Pareto"),
				    type=c("primary","qq"), ...) {
    args <- list(...)

    type <- match.arg(type)

    if(!tform) {
        if(is.element(model, c("Weibull","Frechet"))) mod2 <- "GEV"
        else mod2 <- model
    
        if(!is.element(model,c("PP","GP","Beta","Pareto"))) yq <- qevd(xp, loc=loc, scale=scale, shape=shape, type=mod2)
        else yq <- qevd(xp, threshold=u, loc=loc, scale=scale, shape=shape, type=mod2)
        if(is.null(args$main)) {
            if(type=="primary") m2 <- ""
            else m2 <- deparse(x$call)
            if(is.element(model, c("GEV","Weibull","Frechet","Gumbel"))) plot(yq, sort(y), xlab="Model Quantiles", ylab="Empirical Quantiles", main=m2)
            else plot(yq, sort(y[eid]), xlab="Model Quantiles", ylab="Empirical Quantiles", main=m2, ...)
        } else {
            if(is.element(model, c("GEV","Weibull","Frechet","Gumbel"))) plot(yq, sort(y), xlab="Model Quantiles", ylab="Empirical Quantiles", ...)
            else plot(yq, sort(y[eid]), xlab="Model Quantiles", ylab="Empirical Quantiles", ...)
        }
    } else {
	if(is.null(args$main)) {
                if(type=="primary") m2 <- ""
                else m2 <- deparse(x$call)

                if(is.element(model, c("GEV","Weibull","Frechet"))) m2 <- paste(m2, "(Gumbel Scale)", sep="\n")
                else if(is.element(model, c("PP", "GP", "Beta", "Pareto"))) m2 <-  paste(m2, "Exponential Scale", sep="\n")

                if(is.element(model, c("GEV","Weibull","Gumbel","Frechet"))) plot(-log(-log(sort(xp))), sort(ytrans), main=m2,
                                                                                    xlab="(Standardized) Model Quantiles", ylab="Empirical Residual Quantiles", ...)
                else if(is.element(model, c("GP","Beta","Exponential","Pareto"))) plot(-log(1 - xp), sort(ytrans), main=m2,
                                                                                    xlab="(Standardized) Residual Quantiles", ylab="Empirical Residual Quantiles", ...)
                else if(model=="PP") plot(-log(1 - xp), sort(-log(ytrans)), main=m2, xlab="(Standardized) Residual Quantiles", ylab="Empirical Residual Quantiles", ...)
            } else {
                if(is.element(model, c("GEV","Weibull","Gumbel","Frechet"))) plot(-log(-log(sort(xp))), sort(ytrans), xlab="Model", ylab="Empirical", ...)
                else if(is.element(model, c("GP","Beta","Exponential","Pareto"))) plot(-log(1 - xp), sort(ytrans), xlab="Model", ylab="Empirical", ...)
                else if(model=="PP") plot(-log(1 - xp), sort(-log(ytrans)), xlab="Model", ylab="Empirical", ...)
            }
    }
    abline(0,1)
    invisible()
} # end of 'quantquant.plot.evd' funciton.

quantquant2.plot.evd <- function(x, y, eid, model=c("GEV", "GP", "PP", "Gumbel", "Frechet", "Weibull", "Exponential", "Beta", "Pareto"), type=c("primary","qq2"), ...) {

    args <- list(...)

    type <- match.arg(type)

    if(is.fixedfevd(x)) z <- rextRemes(x, n=x$n)
    else z <- rextRemes(x)

    if(!is.element(model, c("PP","GP","Beta","Exponential","Pareto"))) {

        yQQ <- y
        if(is.null(args$xlab)) xl <- paste(x$data.name[1], " Empirical Quantiles", sep="")
        else xl <- args$xlab

    } else {

        yQQ <- y[eid]
        if(is.null(args$xlab)) {

            if(length(x$threshold)==1) xl <- paste(x$data.name[1], "( > ", x$threshold, ") Empirical Quantiles", sep="")
            else xl <- paste(x$data.name[1], "( > threshold) Empirical Quantiles", sep="")

        } else xl <- args$xlab
    }

    if(!is.null(args$main)) {

        if(type=="primary") mQQ <- ""
        else mQQ <- deparse(x$call)
        qqplot(yQQ, z, main=mQQ, xlab=xl, ylab="Quantiles from Model Simulated Data")

    } else {

        qqplot(yQQ, z, xlab=xl, ylab="Quantiles from Model Simulated Data")

    }

    invisible()

} # end of 'quantquant2.plot.evd' function.

histplot.evd <- function(x, y, u=u, loc=loc, scale=scale, shape=shape, ytrans=NULL, tform=FALSE, eid,
			    model= c("GEV", "GP", "PP", "Gumbel", "Frechet", "Weibull", "Exponential", "Beta", "Pareto"), hist.args=NULL, ...) {
    args <- list(...)
    if(!tform) {
            if(is.null(args$main)) {
                if(model != "PP") m4 <- paste(deparse(x$call), "Histogram", sep="\n")
                else m4 <- paste(deparse(x$call), "\n", paste("Histogram (", x$period.basis, " maxima)", sep=""))
            }
            if(is.element(model, c("GP", "Beta", "Exponential", "Pareto"))) yh <- y[eid]
            else if(model != "PP") yh <- y
            else {
                blocks <- rep(1:x$span, each=x$npy)
                n2 <- length(blocks)
                if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
                else if(n2 > x$n) blocks <- blocks[1:x$n]
                yh <- c(aggregate(y, by=list(blocks), max)$x)
            }
        } else { # Eric 8/14/13
	    if(is.element(model, c("PP", "GP", "Exponential", "Beta", "Pareto"))) stop("plot.fevd: invalid type argument for this model.")
            if(is.null(args$main)) m4 <- paste(deparse(x$call), "Histogram of Transformed Data", sep="\n")
            # if(model != "PP") yh <- ytrans
	    # else {
	# 	blocks <- rep(1:x$span, each=x$npy)
         #        n2 <- length(blocks)
          #       if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
           #      else if(n2 > x$n) blocks <- blocks[1:x$n]
            #     yh <- c(aggregate(ytrans, by=list(blocks), max)$x)
	    # } # end of if else 'model != PP' stmts.
        }

   if(is.null(hist.args)) {
       if(!is.null(args$ylim)) {
           if(is.null(args$main)) {
                    if(is.null(args$col)) hist(yh, col="darkblue", freq=FALSE, breaks="FD", xlab=x$data.name[1], main=m4, ...)
                    else  hist(yh, freq=FALSE, breaks="FD", main=m4, xlab=x$data.name[1], ...)
           } else {
                    if(is.null(args$col)) hist(yh, col="darkblue", freq=FALSE, breaks="FD", xlab=x$data.name[1], ...)
                    else  hist(yh, freq=FALSE, breaks="FD", xlab=x$data.name[1], ...)
                }
           } else {
                if(is.null(args$main)) {
                   if(is.null(args$col)) hist(yh, col="darkblue", freq=FALSE, breaks="FD", xlab=x$data.name[1], main=m4, ylim=c(0,1.5), ...)
                   else  hist(yh, freq=FALSE, breaks="FD", main=m4, xlab=x$data.name[1], ylim=c(0,1.5), ...)
            } else {
               if(is.null(args$col)) hist(yh, col="darkblue", freq=FALSE, breaks="FD", xlab=x$data.name[1], ylim=c(0,1.5), ...)
               else  hist(yh, freq=FALSE, breaks="FD", xlab=x$data.name[1], ylim=c(0,1.5), ...)
            }
        }
    } else do.call("hist", c(list(yh), hist.args))
    xh <- seq(min(yh, na.rm=TRUE), max(yh, na.rm=TRUE),,100)
    if(is.element(model, c("Gumbel", "Weibull", "Frechet"))) mod2 <- "GEV"
    else if(is.element(model, c("PP", "Pareto", "Frechet", "Beta", "Exponential"))) mod2 <- "GP"
    else mod2 <- model

    if(mod2 != "GP") ymod <- devd(xh, loc=loc, scale=scale, shape=shape, threshold=u, type=mod2)
    else ymod <- devd(xh - u, loc=loc, scale=scale, shape=shape, threshold=u, type=mod2)

    lines(xh, ymod, lty=2, col="blue", lwd=1.5)

    invisible()

} # end of 'hitsplot.evd' function.

densplot.evd <- function(x, y, u, loc, scale, shape, tform=FALSE, eid, ytrans=NULL,
			    model= c("GEV", "GP", "PP", "Gumbel", "Frechet", "Weibull", "Exponential", "Beta", "Pareto"), density.args=NULL,
			    type=c("primary","density"), leg=TRUE, ...) {

    args <- list(...)

    type <- match.arg(type)

    if(!tform) {

            if(!is.element(model, c("PP","GP","Beta","Exponential","Pareto"))) yd <- y
            else if(model != "PP") yd <- y[eid]
            else {

                blocks <- rep(1:x$span, each=x$npy)
                n2 <- length(blocks)
                if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
                else if(n2 > x$n) blocks <- blocks[1:x$n]

		ytmp <- c(datagrabber(x, cov.data = FALSE))
		ytmp2 <- c(aggregate(ytmp, by = list(blocks), max)$x)

		mbid <- logical(x$n)

		for(i in 1:(length(unique(blocks)))) {

		    tmpind <- blocks == blocks[i]
		    tmpind2 <- ytmp == ytmp2[i]
		    tmpind2[is.na(tmpind2)] <- FALSE

		    finind <- tmpind & tmpind2

		    if(sum(finind) > 1) {

			numsind <- (1:x$n)[finind]
			numsind <- numsind[1]
			finind[-numsind] <- FALSE

		    }

		    mbid[finind] <- TRUE

		} # end of for 'i' loop.

		yd <- y[mbid]
                # yd <- c(aggregate(y, by=list(blocks), max)$x)
            }

            if(is.null(density.args)) yd <- density(yd)
            else yd <- do.call("density", c(list(yd), density.args))

            if(is.null(args$ylim)) {
                    yld <- range(yd$y, finite=TRUE)
                    yld[1] <- min(yld[1], 0)
                    # yld[2] <- max(yld[2]+0.5, 1)
            }

            if(is.null(args$main)) {
                if(type=="primary") m3 <- ""
                else m3 <- deparse(x$call)
                if(is.null(args$ylim)) plot(yd, main=m3, ylim=yld, ...)
                else plot(yd, main=m3, ...)
            } else {
                if(is.null(args$ylim)) {
                    plot(yd, ylim=yld, ...)
                } else plot(yd, ...)
            }
	    if(is.element(model, c("PP", "Gumbel", "Weibull", "Frechet"))) mod2 <- "GEV"
    	    else if(is.element(model, c("Pareto", "Frechet", "Beta", "Exponential"))) mod2 <- "GP"
    	    else mod2 <- model

    	    xd <- seq(min(yd$x, na.rm=TRUE), max(yd$x, na.rm=TRUE),,100)
    	    if(mod2 != "GP") ymod <- devd(xd, loc=loc, scale=scale, shape=shape, threshold=u, type=mod2)
	    else ymod <- devd(xd - u, loc=loc, scale=scale, shape=shape, threshold=u, type=mod2)

    	    lines(xd, ymod, lty=2, col="blue", lwd=1.5)
	    if(leg) legend("topright", legend=c("Data","Model"), col=c("black","blue"), lty=c(1,2), lwd=c(1,1.5), bty="n")

        } else {

	    # if(model != "PP") ytrans2 <- ytrans
	    # else {

	# 	blocks <- rep(1:x$span, each=x$npy)
         #        n2 <- length(blocks)
          #       if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
           #      else if(n2 > x$n) blocks <- blocks[1:x$n]
            #     ytrans2 <- c(aggregate(ytrans, by=list(blocks), max)$x)

	    # }

            if(is.null(density.args)) yd <- density(ytrans)
            else yd <- do.call("density", c(list(ytrans), density.args))

            if(is.null(args$ylim)) {
                    yld <- range(yd$y, finite=TRUE)
                    yld[1] <- min(yld[1], 0)
                    # yld[2] <- max(yld[2]+0.5, 1)
            }

            if(is.null(args$main)) {
                if(type=="primary") m3 <- ""
                else m3 <- paste(deparse(x$call), "Transformed", sep="\n")
                if(is.null(args$ylim)) plot(yd, main=m3, ylim=yld, ...)
                else plot(yd, main=m3, ...)
            } else {
                if(is.null(args$ylim)) plot(yd, ylim=yld, ...)
                else plot(yd, ...)
            }
	    if(is.element(model, c("Gumbel", "Weibull", "Frechet"))) mod2 <- "GEV"
    	    else if(is.element(model, c("PP", "Pareto", "Frechet", "Beta", "Exponential"))) mod2 <- "GP"
    	    else mod2 <- model

    	    xd <- seq(min(yd$x, na.rm=TRUE), max(yd$x, na.rm=TRUE),,100)
    	    ymod <- devd(xd, loc=0, scale=1, shape=0, threshold=0, type=mod2)
    	    lines(xd, ymod, lty=2, col="blue", lwd=1.5)

	    if(leg) legend("topright", legend=c("Transformed Data","Standardized Model"), col=c("black","blue"), lty=c(1,2), lwd=c(1,1.5), bty="n")
        } # end of if else '!tform' stmts.
        # if(!is.element(model, c("pp","gp","beta","exponential","pareto"))) points(y, rep(0, length(y)), pch="|", cex=0.5)
        # else points(y[eid], rep(0, length(y[eid])), pch="|", cex=0.5)

    invisible()

} # end of 'densplot.evd' function.

rlplot.evd <- function(x, xp, y, u, eid, rperiods, tform=FALSE, model= c("GEV", "GP", "PP", "Gumbel", "Frechet", "Weibull", "Exponential", "Beta", "Pareto"),
			    type=c("primary","rl"), leg=TRUE, a=0, ...) {
    args <- list(...)

    type <- match.arg(type)

    if(x$method != "Lmoments") { # Eric -- 8/27/13

        const.thresh <- check.constant(x$par.models$threshold)
        const.loc <- check.constant(x$par.models$location)
        const.scale <- check.constant(x$par.models$scale)
        const.shape <- check.constant(x$par.models$shape)

        if(is.element(model, c("PP", "GP", "Exponential", "Beta", "Pareto")) && !const.thresh && all(c(const.loc, const.scale, const.shape)) && type == "rl")
            stop("rlplot.evd: invalid type argument for POT models with varying thresholds but constant parameters (are you sure you about this model choice?).")

    }

    if(is.null(args$main)) {
        if(type=="primary") m5 <- ""
        else m5 <- deparse(x$call)
    }
    if(model=="PP") {

        if(!tform) mod2 <- "GEV"
	else mod2 <- "GEV"

        if(is.null(args$main)) {
	    if(!tform) {
                if(type=="rl") m5 <- paste(m5, "Return Levels based on approx. equivalent GEV", sep="\n")
                else m5 <- paste("Return Levels based on approx.", "equivalent GEV", sep="\n")
	    } else {
		if(type=="rl") m5 <- paste(m5, "Return Levels based on approx. equivalent GEV df", sep="\n")
		else m5 <- paste("Return Levels based on approx.", "equivalent GEV df", sep="\n")
	    }
        }

        blocks <- rep(1:x$span, each=x$npy)
        n2 <- length(blocks)
        if(n2 < x$n) blocks <- c(blocks, rep(blocks[n2], x$n - n2))
        else if(n2 > x$n) blocks <- blocks[1:x$n]
        yEmp <- c(aggregate(y, by=list(blocks), max)$x)

    } else mod2 <- model

    if(!tform) {

        bds <- ci(x, return.period=rperiods)
	yrl <- bds[,2]
        if(is.null(args$ylim)) yl <- range(c(bds), finite=TRUE)
            
        if(is.element(model, c("PP", "GEV", "Gumbel", "Weibull", "Frechet"))) xrl <- -1/(log(1 - 1/rperiods))
        else xrl <- rperiods

        if(is.null(x$units)) ylb <- "Return Level"
        else ylb <- paste("Return Level (", x$units, ")", sep="")

        xlb <- paste("Return Period (", x$period.basis, "s)", sep="")

        if(is.null(args$main)) {
            if(!is.null(args$ylim)) plot(xrl, yrl, type="l", log="x", xlab=xlb, ylab=ylb, main=m5, ...)
            else plot(xrl, yrl, type="l", log="x", ylim=yl, xlab=xlb, ylab=ylb, main=m5, ...)
        } else {
            if(!is.null(args$ylim)) plot(xrl, yrl, type="l", log="x", xlab=xlb, ylab=ylb, ...)
            else plot(xrl, yrl, type="l", log="x", ylim=yl, xlab=xlb, ylab=ylb, ...)
        }

        lines(xrl, bds[,1], col="gray", lty=2, lwd=2)
        lines(xrl, bds[,3], col="gray", lty=2, lwd=2)

        if(is.element(model, c("GEV", "Gumbel", "Weibull", "Frechet"))) points(-1/log(xp), sort(y))
        else if(is.element(model, c("GP", "Beta", "Pareto", "Exponential"))) {

            n2 <- x$n
            if(is.null(a)) xp2 <- ppoints(n2)
            else xp2 <- ppoints(n2, a=a)
            sdat <- sort(y)
            points(-1/log(xp2)[sdat > u]/x$npy, sdat[sdat > u])

        } else if(model == "PP") {

            if(is.null(a)) xp2 <- ppoints(length(yEmp))
            else xp2 <- ppoints(length(yEmp), a=a)
            points(-1/log(xp2), sort(yEmp))

        }

    } else {

	np <- length(rperiods)
	n <- length(y)
        effrl <- matrix(NA, n, np)
        for(i in 1:np) effrl[,i] <- return.level(x, return.period=rperiods[i])

        if(is.null(args$ylim)) {

                yl <- range(c(y, c(effrl)), finite=TRUE)
                yl[2] <- yl[2] + sign(yl[2]) * log2(abs(yl[2]))

                if(is.null(args$main)) {

                    if(!is.element(model, c("GP", "Beta", "Exponential", "Pareto"))) plot(y, type="l", xlab="index", main=m5, ylim=yl, ...)
		    else plot(y[eid], type = "l", xlab = "index", ylim = yl, main = m5, ...)


        	} else {

                    if(!is.element(model, c("GP", "Beta", "Exponential", "Pareto"))) plot(y, type="l", xlab="index", ...)
                    else plot(y[eid], type="l", xlab="index", ylim=yl, ...)

                }

        } else {

            if(is.null(args$main)) {

                    if(!is.element(model, c("GP", "Beta", "Exponential", "Pareto"))) plot(y, type="l", xlab="index", main=m5, ...)
                    else plot(y[eid], type = "l", xlab = "index", main = m5, ...)

                } else {

                    if(!is.element(model, c("GP", "Beta","Exponential","Pareto"))) plot(y, type="l", xlab="index", ...)
                    else plot(y[eid], type="l", xlab="index", ...)

                }

        } # end of if else 'ylim passed via '...' stmts.

        for(i in 1:np) lines(effrl[,i], lty=i, col=i+1)

	if(is.element(model, c("PP", "GP", "Beta", "Exponential", "Pareto"))) {

	    if(length(u) > 1 && all(u == u[1])) u <- u[1]
	    if(length(u) == 1) abline(h=u, col="darkorange")
	    else lines(u, col="darkorange")
	    if(leg) legend("topleft", legend=c(paste(rperiods, "-", x$period.basis, " level", sep=""), "threshold"), lty=c(1:np,1), col=c(2:(np+1),"darkorange"), bg="white")

	} else if(leg) legend("topleft", legend=paste(rperiods, "-", x$period.basis, " level", sep=""), lty=1:np, col=2:(np+1), bg="white")

    } # end of if else '!tform' stmts.

    invisible()

} # end of 'rlplot.evd' function.

eeplot <- function(x, type = c("both", "Zplot", "Wplot"), set.pw = FALSE, d = NULL, ...) {

    type <- match.arg(type)

    theCall <- match.call()

    if(x$type != "PP") stop("eeplot: Sorry, this funciton is for PP fits only.")

    y <- datagrabber(x, cov.data = FALSE)
    u <- x$threshold
    p <- findpars(x)

    Tk <- (1:x$n)[y > u]

    out <- Tk

    if(set.pw) {

        if(type == "both") par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
        else par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))

    }

    if(is.element(type, c("both", "Zplot"))) {

        lambda <- (1 + p$shape * (u - p$location)/p$scale)^(-1/p$shape)

        if(is.null(d)) {
            
            if(x$time.units == "days") {

                if(x$period.basis == "year") d <- 365.25
                else stop("eeplot: invalid period basis, perhaps set the d argument.")

            } else if(x$time.units == "months") {

                if(x$period.basis == "year") d <- 12
                else stop("eeplot: invalid period basis, perhaps set the d argument.")

            } else if(x$time.units == "years") {

                if(x$period.basis == "year") d <- 1
                else stop("eeplot: invalid period basis, perhaps set the d argument.")

            } else if(x$time.units == "hours") {

                if(x$period.basis == "year") d <- 8766 # 24 * 365.25
                else stop("eeplot: invalid period basis, perhaps set the d argument.")

            } else if(x$time.units == "minutes") {

                if(x$period.basis == "year") d <- 525960 # 60 * 24 * 365.25
                else stop("eeplot: invalid period basis, perhaps set the d argument.")

            } else if(x$time.units == "seconds") {

                if(x$period.basis == "year") d <- 31557600 # 60 * 60 * 24 * 365.25
                else stop("eeplot: invalid period basis, perhaps set the d argument.")

            } else {

		tmp.units <- unlist(strsplit(x$time.units, split="/"))
                if(length(tmp.units) != 2) stop("fevd: invalid time.units argument.")
                numper <- as.numeric(tmp.units[1])
                if(is.na(numper)) stop("fevd: invalid time.units argument.")
                pertiid <- tmp.units[2]

                if(!is.element(pertiid, c("day","month","year","hour","minute","second"))) stop("fevd: invalid time.units argument.")

                if(pertiid=="year") d <- numper
                else if(pertiid=="month") d <- numper*12
                else if(pertiid=="day") d <- numper*365.25
		else if(pertiid == "hour") d <- numper * 24 * 365.25
		else if(pertiid == "minute") d <- numper * 60 * 24 * 365.25
		else if(pertiid == "second") d <- numper * 60 * 60 * 24 * 365.25

	    }

        }

        lambda <- lambda/d
	if(length(lambda) == 1) lambda <- rep(lambda, x$n)

        Zk <- diff(c(0, cumsum(lambda)[Tk]))
	out <- cbind(out, Zk)
	colnames(out) <- c("Exceed Time", "Zk")

        qqPlot(Zk, distribution = "exp", xlab = "Expected Values\nUnder exponential(1)", ylab = "Observed Z_k Values",
		col.lines = "grey", grid = FALSE, ...)

	abline(0, 1, col = "darkorange", lty = 2)

	legend("topleft", legend = c("1-1 line", "regression line", "95% confidence bands"), col = c("darkorange", "grey", "grey"), lty = c(2, 1, 2))

    }

    if(is.element(type, c("both", "Wplot"))) {

        if(length(u) == 1) u <- rep(u, x$n)

        Wk <- (1/p$shape[Tk]) * log(1 + p$shape[Tk] * (y[Tk] - u[Tk])/(p$scale[Tk] + p$shape[Tk] * (u[Tk] - p$location[Tk])))

        qqPlot(Wk, distribution = "exp", xlab = "Expected Values\nUnder exponential(1)", ylab = "Observed W_k Values", 
		col.lines = "grey", grid = FALSE, ...)

	abline(0, 1, col = "darkorange", lty = 2)

        legend("topleft", legend = c("1-1 line", "regression line", "95% confidence bands"), col = c("darkorange", "grey", "grey"), lty = c(2, 3, 1))

	out <- cbind(out, Wk)
	if(is.matrix(out)) colnames(out) <- c("Exceed Time", "Zk", "Wk")
	else colnames(out) <- c("Exceed Time", "Wk")
    }

    if(set.pw) mtext(deparse(x$call), line=0.5, outer=TRUE)

    attr(out, "call") <- theCall
    class(out) <- "ee"
    invisible(out)

} # end of 'eeplot' function.

plot.ee <- function(x, pch2 = "+", col2 = "gray", xlab = "Time of Exceeding Threshold", ylab = "", ...) {

    args <- list(...)

    if(is.null(args$pch)) pch1 <- "o"
    else pch1 <- args$pch
            
    if(is.null(args$col)) col1 <- 1
    else col1 <- args$col


    if(is.null(args$ylim)) {

	yl <- range(x[,-1], finite = TRUE)
	plot(x[,1], x[,2], ylim = yl, xlab = xlab, ylab = ylab, ...)

    } else plot(x[,1], x[,2], xlab = xlab, ylab = ylab, ...)

    if(dim(x)[2] == 3) {

        points(x[,1], x[,3], pch = pch2, col = col2)
        legend("topleft", legend = c("Zk", "Wk"), pch = c(pch1, pch2), col = c(col1, col2), bty = "n")

    }

    if(is.null(args$main)) title(attributes(x)$call)

    invisible()

} # end of 'plot.ee' function.
