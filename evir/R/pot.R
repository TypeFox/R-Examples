"gpd" <- 
function(data, threshold = NA, nextremes = NA, method = c("ml","pwm"),
         information = c("observed","expected"), ...)
{
    data <- as.numeric(data)
    n <- length(data)
    if(is.na(nextremes) && is.na(threshold))
        stop("Enter either a threshold or the number of upper extremes")
    if(!is.na(nextremes) && !is.na(threshold))
        stop("Enter EITHER a threshold or the number of upper extremes")
    if(!is.na(nextremes))
        threshold <- findthresh(data, nextremes)
    exceedances <- data[data > threshold]
    excess <- exceedances - threshold
    Nu <- length(excess)
    xbar <- mean(excess)
    method <- match.arg(method)
    if(method == "ml") {
        s2 <- var(excess)
        xi0 <- -0.5 * (((xbar * xbar)/s2) - 1)
        beta0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
        theta <- c(xi0, beta0)
        negloglik <- function(theta, tmp)
        {
       	    xi <- theta[1]
            beta <- theta[2]
	    cond1 <- beta <= 0
	    cond2 <- (xi <= 0) && (max(tmp) > ( - beta/xi))
	    if(cond1 || cond2)
	  	f <- 1e+06
	    else {
	    	y <- logb(1 + (xi * tmp)/beta)
	        y <- y/xi
	        f <- length(tmp) * logb(beta) + (1 + xi) * sum(y)
	    }
	    f
	}
        fit <- optim(theta, negloglik, hessian = TRUE, ..., tmp = excess)
        if(fit$convergence)
            warning("optimization may not have succeeded")
        par.ests <- fit$par
        converged <- fit$convergence
        nllh.final <- fit$value
        information <- match.arg(information)
        if(information == "observed") varcov <- solve(fit$hessian)
        if(information == "expected") {
            one <- (1 + par.ests[1])^2 / Nu
	    two <- (2 * (1 + par.ests[1]) * par.ests[2]^2) / Nu
	    cov <-  - ((1 + par.ests[1]) * par.ests[2]) / Nu
	    varcov <- matrix(c(one, cov, cov, two), 2)
	}
    }
    if(method == "pwm") {
        a0 <- xbar
        gamma <- -0.35
        delta <- 0
        pvec <- ((1:Nu) + gamma)/(Nu + delta)
        a1 <- mean(sort(excess) * (1 - pvec))
        xi <- 2 - a0/(a0 - 2 * a1)
        beta <- (2 * a0 * a1)/(a0 - 2 * a1)
        par.ests <- c(xi, beta)
        denom <- Nu * (1 - 2 * xi) * (3 - 2 * xi)
        if(xi > 0.5) {
            denom <- NA
      	    warning("Asymptotic standard errors not available for",
                    "PWM Method when xi > 0.5")
        }
        one <- (1 - xi) * (1 - xi + 2 * xi^2) * (2 - xi)^2
        two <- (7 - 18 * xi + 11 * xi^2 - 2 * xi^3) * beta^2
        cov <- beta * (2 - xi) * (2 - 6 * xi + 7 * xi^2 - 2 * xi^3)
        varcov <- matrix(c(one, cov, cov, two), 2) / denom
        information <- "expected"
        converged <- NA
        nllh.final <- NA
    }
    par.ses <- sqrt(diag(varcov))
    p.less.thresh <- 1 - Nu/n
    out <- list(n = length(data), data = exceedances, threshold =
        threshold, p.less.thresh = p.less.thresh, n.exceed = Nu,
        method = method, par.ests = par.ests, par.ses = par.ses,
        varcov = varcov, information = information, converged =
        converged, nllh.final = nllh.final)
    names(out$par.ests) <- c("xi", "beta")
    names(out$par.ses) <- c("xi", "beta")
    class(out) <- "gpd"
    out
}

"gpd.q" <- 
function(x, pp, ci.type = c("likelihood","wald"), ci.p = 0.95,
         like.num = 50)
{
    if(x$dist != "gpd")
        stop("This function is used only with GPD curves")
    if(length(pp) > 1)
	stop("One probability at a time please")
    threshold <- x$lastfit$threshold
    par.ests <- x$lastfit$par.ests
    xihat <- par.ests["xi"]
    betahat <- par.ests["beta"]
    varcov <- x$lastfit$varcov
    p.less.thresh <- x$lastfit$p.less.thresh
    lambda <- 1
    if(x$type == "tail") lambda <- 1/(1 - p.less.thresh)
    a <- lambda * (1 - pp)
    gfunc <- function(a, xihat) (a^( - xihat) - 1) / xihat
    gfunc.deriv <- function(a, xihat)
        ( - (a^( - xihat) - 1)/xihat - a^( - xihat) * logb(a)) / xihat
    q <- threshold + betahat * gfunc(a, xihat)
    if(q < x$plotmax) abline(v = q, lty = 2)
    out <- as.numeric(q)
    ci.type <- match.arg(ci.type)
    if(ci.type == "wald") {
        if(class(x$lastfit) != "gpd")
	    stop("Wald method requires model be fitted with gpd (not pot)")
	scaling <- threshold
	betahat <- betahat/scaling
	xivar <- varcov[1, 1]
        betavar <- varcov[2, 2]/(scaling^2)
	covar <- varcov[1, 2]/scaling
	term1 <- betavar * (gfunc(a, xihat))^2
	term2 <- xivar * (betahat^2) * (gfunc.deriv(a, xihat))^2
	term3 <- 2 * covar * betavar * gfunc(a, xihat) * gfunc.deriv(a, xihat)
	qvar <- term1 + term2 + term3
	if(qvar < 0) stop("Negative estimate of quantile variance")
	qse <- scaling * sqrt(qvar)
	qq <- qnorm(1 - (1 - ci.p)/2)
	upper <- q + qse * qq
	lower <- q - qse * qq
        if(upper < x$plotmax) abline(v = upper, lty = 2, col = 2)
	if(lower < x$plotmax) abline(v = lower, lty = 2, col = 2)
	out <- as.numeric(c(lower, q, qse, upper))
	names(out) <- c("Lower CI", "Estimate", "Std.Err", "Upper CI")
    }
    if(ci.type == "likelihood") {
        parloglik <- function(theta, tmp, a, threshold, xpi)
	{
	    beta <- (theta * (xpi - threshold))/(a^( - theta) - 1)
	    if((beta <= 0) || ((theta <= 0) && (max(tmp) > ( - beta/theta))))
	        f <- 1e+06
	    else {
	        y <- logb(1 + (theta * tmp)/beta)
		y <- y/theta
		f <- length(tmp) * logb(beta) + (1 + theta) * sum(y)
	    }
	    f
	}
	theta <- xihat
	parmax <- NULL
	xp <- exp(seq(from = logb(threshold), to = logb(x$plotmax),
                      length = like.num))
        excess <- as.numeric(x$lastfit$data - threshold)
	for(i in 1:length(xp)) {
            fit2 <- optim(theta, parloglik, method = "BFGS", hessian = FALSE,
                tmp = excess, a = a, threshold = threshold, xpi = xp[i])
	    parmax <- rbind(parmax, fit2$value)
	}
	parmax <-  - parmax
	overallmax <-  - parloglik(xihat, excess, a, threshold, q)
	crit <- overallmax - qchisq(0.999, 1)/2
	cond <- parmax > crit
	xp <- xp[cond]
	parmax <- parmax[cond]
	par(new = TRUE)
	dolog <- ""
        if(x$alog == "xy" || x$alog == "x") dolog <- "x"
	plot(xp, parmax, type = "n", xlab = "", ylab = "", axes = FALSE,
	     xlim = range(x$plotmin, x$plotmax),
	     ylim = range(overallmax, crit), log = dolog)
	axis(4, at = overallmax - qchisq(c(0.95, 0.99), 1)/2,
             labels = c("95", "99"), tick = TRUE)
	aalpha <- qchisq(ci.p, 1)
	abline(h = overallmax - aalpha/2, lty = 2, col = 2)
	cond <- !is.na(xp) & !is.na(parmax)
	smth <- spline(xp[cond], parmax[cond], n = 200)
	lines(smth, lty = 2, col = 2)
	ci <- smth$x[smth$y > overallmax - aalpha/2]
	out <- c(min(ci), q, max(ci))
	names(out) <- c("Lower CI", "Estimate", "Upper CI")
    }
    out
}

"gpd.sfall" <- 
function(x, pp, ci.p = 0.95, like.num = 50)
{
    if(x$dist != "gpd")
       	stop("This function is used only with GPD curves")
    if(length(pp) > 1)
       	stop("One probability at a time please")
    threshold <- x$lastfit$threshold
    par.ests <- x$lastfit$par.ests
    xihat <- par.ests["xi"]
    betahat <- par.ests["beta"]
    varcov <- x$lastfit$varcov
    p.less.thresh <- x$lastfit$p.less.thresh
    lambda <- 1
    if(x$type == "tail") lambda <- 1/(1 - p.less.thresh)
    a <- lambda * (1 - pp)
    gfunc <- function(a, xihat) (a^( - xihat) - 1) / xihat
    q <- threshold + betahat * gfunc(a, xihat)
    s <- q + (betahat + xihat * (q - threshold))/(1 - xihat)
    if(s < x$plotmax) abline(v = s, lty = 2)
    out <- as.numeric(s)
    parloglik <- function(theta, tmp, a, threshold, xpi)
    {
        beta <- ((1 - theta) * (xpi - threshold)) /
          (((a^( - theta) - 1)/theta) + 1)
	if((beta <= 0) || ((theta <= 0) && (max(tmp) > ( - beta/theta))))
	    f <- 1e+06
	else {
	    y <- logb(1 + (theta * tmp)/beta)
	    y <- y/theta
	    f <- length(tmp) * logb(beta) + (1 + theta) * sum(y)
	}
	f
    }
    theta <- xihat
    parmax <- NULL
    xp <- exp(seq(from = logb(threshold), to = logb(x$plotmax),
                  length = like.num))
    excess <- as.numeric(x$lastfit$data - threshold)
    for(i in 1:length(xp)) {
        fit2 <- optim(theta, parloglik, method = "BFGS", hessian = FALSE,
                tmp = excess, a = a, threshold = threshold, xpi = xp[i])
        parmax <- rbind(parmax, fit2$value)
    }
    parmax <-  - parmax
    overallmax <-  - parloglik(xihat, excess, a, threshold, s)
    crit <- overallmax - qchisq(0.999, 1)/2
    cond <- parmax > crit
    xp <- xp[cond]
    parmax <- parmax[cond]
    par(new = TRUE)
    dolog <- ""
    if(x$alog == "xy" || x$alog == "x") dolog <- "x"
    plot(xp, parmax, type = "n", xlab = "", ylab = "", axes = FALSE,
         xlim = range(x$plotmin, x$plotmax), ylim =
         range(overallmax, crit), log = dolog)
    axis(4, at = overallmax - qchisq(c(0.95, 0.99), 1)/2,
         labels = c("95", "99"), tick = TRUE)
    aalpha <- qchisq(ci.p, 1)
    abline(h = overallmax - aalpha/2, lty = 2, col = 2)
    cond <- !is.na(xp) & !is.na(parmax)
    smth <- spline(xp[cond], parmax[cond], n = 200)
    lines(smth, lty = 2, col = 2)
    ci <- smth$x[smth$y > overallmax - aalpha/2]
    out <- c(min(ci), s, max(ci))
    names(out) <- c("Lower CI", "Estimate", "Upper CI")
    out
}

"plot.gpd" <- 
function(x, optlog = NA, extend = 1.5, labels = TRUE, ...)
{
    data <- as.numeric(x$data)
    threshold <- x$threshold
    xi <- x$par.ests["xi"]
    beta <- x$par.ests["beta"]
    choices <- c("Excess Distribution", "Tail of Underlying Distribution",
      	"Scatterplot of Residuals", "QQplot of Residuals")
    tmenu <- paste("plot:", choices)
    pick <- 1
    lastcurve <- NULL
    while(pick > 0) {
        pick <- menu(tmenu, title =
                     "\nMake a plot selection (or 0 to exit):")
        if(pick >= 3) {
            excess <- data - threshold
       	    res <- logb(1 + (xi * excess)/beta) / xi
            lastcurve <- NULL
	}
        if(pick == 3) {
      	    plot(res, ylab = "Residuals", xlab = "Ordering", ...)
	    lines(lowess(1:length(res), res))
	}
        if(pick == 4) qplot(res, ...)
        if(pick == 1 || pick == 2) {
            plotmin <- threshold
     	    if(extend <= 1) stop("extend must be > 1")
	    plotmax <- max(data) * extend
            xx <- seq(from = 0, to = 1, length = 1000)
      	    z <- qgpd(xx, xi, threshold, beta)
       	    z <- pmax(pmin(z, plotmax), plotmin)
       	    ypoints <- ppoints(sort(data))
       	    y <- pgpd(z, xi, threshold, beta)
	}
	if(pick == 1) {
	    type <- "eplot"
	    if(!is.na(optlog))
                alog <- optlog
       	    else alog <- "x"
	    if(alog == "xy")
	        stop("Double log plot of Fu(x-u) does\nnot make much sense")
	    yylab <- "Fu(x-u)"
       	    shape <- xi
	    scale <- beta
	    location <- threshold
	}
	if(pick == 2) {
	    type <- "tail"
	    if(!is.na(optlog))
	        alog <- optlog
	    else alog <- "xy"
	    prob <- x$p.less.thresh
	    ypoints <- (1 - prob) * (1 - ypoints)
	    y <- (1 - prob) * (1 - y)
	    yylab <- "1-F(x)"
	    shape <- xi
	    scale <- beta * (1 - prob)^xi
	    location <- threshold - (scale * ((1 - prob)^( - xi) - 1))/xi
	}
	if(pick == 1 || pick == 2) {
	    plot(sort(data), ypoints, xlim = range(plotmin, plotmax),
                 ylim = range(ypoints, y, na.rm = TRUE), xlab = "",
                 ylab = "", log = alog, axes = TRUE, ...)
	    lines(z[y >= 0], y[y >= 0])
	    if(labels) {
	        xxlab <- "x"
		if(alog == "x" || alog == "xy" || alog == "yx")
		    xxlab <- paste(xxlab, "(on log scale)")
	        if(alog == "xy" || alog == "yx" || alog == "y")
		    yylab <- paste(yylab, "(on log scale)")
		title(xlab = xxlab, ylab = yylab)
	    }
	    details <- paste("threshold = ", format(signif(threshold, 3)),
                             "   xi = ", format(signif(shape, 3)),
                             "   scale = ", format(signif(scale, 3)),
                             "   location = ", format(signif(location, 3)),
                             sep = "")
	    print(details)
	    lastcurve <- list(lastfit = x, type = type, dist = "gpd",
                plotmin = plotmin, plotmax = plotmax, alog = alog,
                location = as.numeric(location), shape = as.numeric(shape),
                scale = as.numeric(scale))
	}
    }
    invisible(lastcurve)
}

"quant" <- 
function(data, p = 0.99, models = 30, start = 15, end = 500,
	reverse = TRUE, ci = 0.95, auto.scale = TRUE, labels = TRUE,
	...)
{
    data <- as.numeric(data)
    n <- length(data)
    if(ci) qq <- qnorm(1 - (1 - ci)/2)
    exceed <- trunc(seq(from = min(end, n), to = start, length = models))
    if(p < 1 - min(exceed)/n) {
        cat("Graph may look strange !! \n\n")
	cat(paste("Suggestion 1: Increase `p' above",
            format(signif(1 - min(exceed)/n, 5)), "\n"))
	cat(paste("Suggestion 2: Increase `start' above ",
            ceiling(length(data) * (1 - p)), "\n"))
    }
    gpd.dummy <- function(nex, data)
    {
	out <- gpd(data = data, nextremes = nex, information = "expected")
	c(out$threshold, out$par.ests[1], out$par.ests[2],
          out$varcov[1, 1], out$varcov[2, 2], out$varcov[1, 2])
    }
    mat <- apply(as.matrix(exceed), 1, gpd.dummy, data = data)
    thresh <- mat[1,  ]
    xihat <- mat[2,  ]
    betahat <- mat[3,  ]
    lambda <- length(data)/exceed
    a <- lambda * (1 - p)
    gfunc <- function(a, xihat) (a^( - xihat) - 1) / xihat
    qest <- thresh + betahat * gfunc(a, xihat)
    l <- u <- qest
    yrange <- range(qest)
    if(ci) {
        xivar <- mat[4,  ]
	betavar <- mat[5,  ]
	covar <- mat[6,  ]
	scaling <- thresh
	betahat <- betahat/scaling
	betavar <- betavar/(scaling^2)
	covar <- covar/scaling
	gfunc.deriv <- function(a, xihat)
	    ( - (a^( - xihat) - 1)/xihat - a^( - xihat) * logb(a)) / xihat
	term1 <- betavar * (gfunc(a, xihat))^2
	term2 <- xivar * (betahat^2) * (gfunc.deriv(a, xihat))^2
	term3 <- 2 * covar * betavar * gfunc(a, xihat) * gfunc.deriv(a, xihat)
	qvar <- term1 + term2 + term3
	if(min(qvar) < 0)
	    stop(paste("Conditioning problems lead to estimated negative",
                       "quantile variance", sep = "\n"))
	qse <- scaling * sqrt(qvar)
	u <- qest + qse * qq
	l <- qest - qse * qq
	yrange <- range(qest, u, l)
    }
    mat <- rbind(thresh, qest, exceed, l, u)
    dimnames(mat) <- list(c("threshold", "qest", "exceedances", "lower",
        "upper"), NULL)
    index <- exceed
    if(reverse) index <-  - exceed
    if(auto.scale)
        plot(index, qest, ylim = yrange, type = "l", xlab = "", ylab = "",
             axes = FALSE, ...)
    else plot(index, qest, type = "l", xlab = "", ylab = "",
              axes = FALSE, ...)
    axis(1, at = index, labels = paste(exceed))
    axis(2)
    axis(3, at = index, labels = paste(format(signif(thresh, 3))))
    box()
    if(ci) {
       	lines(index, l, lty = 2, col = 2)
       	lines(index, u, lty = 2, col = 2)
    }
    if(labels) {
       	labely <- paste(p, "Quantile")
       	if(ci) labely <- paste(labely, " (CI, p = ", ci, ")", sep = "")
	title(xlab = "Exceedances", ylab = labely)
	mtext("Threshold", side = 3, line = 3)
    }
    invisible(mat)
}

"riskmeasures" <- 
function(x, p)
{
    u <- x$threshold
    par.ests <- x$par.ests
    xihat <- par.ests["xi"]
    betahat <- par.ests["beta"]
    p.less.thresh <- x$p.less.thresh
    lambda <- 1/(1 - p.less.thresh)
    quant <- function(pp, xi, beta, u, lambda)
    {
     	a <- lambda * (1 - pp)
       	u + (beta * (a^( - xi) - 1))/xi
    }
    short <- function(pp, xi, beta, u, lambda)
    {
      	a <- lambda * (1 - pp)
       	q <- u + (beta * (a^( - xi) - 1))/xi
       	(q * (1 + (beta - xi * u)/q)) / (1 - xi)
    }
    q <- quant(p, xihat, betahat, u, lambda)
    es <- short(p, xihat, betahat, u, lambda)
    rtn <- cbind(p, quantile = q, sfall = es)
    row.names(rtn) <- NULL
    rtn
}

"shape" <- 
function(data, models = 30, start = 15, end = 500, reverse = TRUE, ci = 
	0.95, auto.scale = TRUE, labels = TRUE, ...)
{
    data <- as.numeric(data)
    qq <- 0
    if(ci) qq <- qnorm(1 - (1 - ci)/2)
    x <- trunc(seq(from = min(end, length(data)), to = start, length = models))
    gpd.dummy <- function(nex, data)
    {
        out <- gpd(data = data, nextremes = nex, information = "expected")
	c(out$threshold, out$par.ests[1], out$par.ses[1])
    }
    mat <- apply(as.matrix(x), 1, gpd.dummy, data = data)
    mat <- rbind(mat, x)
    dimnames(mat) <- list(c("threshold", "shape", "se", "exceedances"), NULL)
    thresh <- mat[1,  ]
    y <- mat[2,  ]
    yrange <- range(y)
    if(ci) {
        u <- y + mat[3,  ] * qq
	l <- y - mat[3,  ] * qq
	yrange <- range(y, u, l)
    }
    index <- x
    if(reverse) index <-  - x
    if(auto.scale)
        plot(index, y, ylim = yrange, type = "l", xlab = "", ylab = "",
	     axes = FALSE, ...)
    else plot(index, y, type = "l", xlab = "", ylab = "", axes = FALSE, ...)
    axis(1, at = index, labels = paste(x), tick = FALSE)
    axis(2)
    axis(3, at = index, labels = paste(format(signif(thresh, 3))), tick = FALSE)
    box()
    if(ci) {
        lines(index, u, lty = 2, col = 2)
	lines(index, l, lty = 2, col = 2)
    }
    if(labels) {
        labely <- "Shape (xi)"
	if(ci) labely <- paste(labely, " (CI, p = ", ci, ")", sep = "")
	title(xlab = "Exceedances", ylab = labely)
	mtext("Threshold", side = 3, line = 3)
    }
    invisible(mat)
}

"tailplot" <- 
function(x, optlog = NA, extend = 1.5, labels = TRUE, ...)
{
    data <- as.numeric(x$data)
    threshold <- x$threshold
    xi <- x$par.ests["xi"]
    beta <- x$par.ests["beta"]
    plotmin <- threshold
    if(extend <= 1) stop("extend must be > 1")
    plotmax <- max(data) * extend
    xx <- seq(from = 0, to = 1, length = 1000)
    z <- qgpd(xx, xi, threshold, beta)
    z <- pmax(pmin(z, plotmax), plotmin)
    ypoints <- ppoints(sort(data))
    y <- pgpd(z, xi, threshold, beta)
    type <- "tail"
    if(!is.na(optlog))
	    alog <- optlog
    else alog <- "xy"
    prob <- x$p.less.thresh
    ypoints <- (1 - prob) * (1 - ypoints)
    y <- (1 - prob) * (1 - y)
    yylab <- "1-F(x)"
    shape <- xi
    scale <- beta * (1 - prob)^xi
    location <- threshold - (scale * ((1 - prob)^( - xi) - 1))/xi
    plot(sort(data), ypoints, xlim = range(plotmin, plotmax), ylim =
         range(ypoints, y, na.rm = TRUE), xlab = "", ylab = "", log = alog,
         axes = TRUE, ...)
    lines(z[y >= 0], y[y >= 0])
    if(labels) {
        xxlab <- "x"
        if(alog == "x" || alog == "xy" || alog == "yx")
	    xxlab <- paste(xxlab, "(on log scale)")
        if(alog == "xy" || alog == "yx" || alog == "y")
	    yylab <- paste(yylab, "(on log scale)")
        title(xlab = xxlab, ylab = yylab)
    }
    lastcurve <- list(lastfit = x, type = type, dist = "gpd",
        plotmin = plotmin, plotmax = plotmax, alog = alog, location = 
	as.numeric(location), shape = as.numeric(shape), scale = 
	as.numeric(scale))
    invisible(lastcurve)
}
