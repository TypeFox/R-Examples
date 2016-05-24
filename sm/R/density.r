"sm.density" <- function(x, h,  model = "none", weights = NA, group = NA, ...) {

    x.name <- deparse(substitute(x))

    data   <- sm.check.data(x, NA, weights, group, ...)
    x      <- data$x
    weights<- data$weights
    group  <- data$group
    nobs   <- data$nobs
    ndim   <- data$ndim
    opt    <- data$options
    
    if (ndim > 3) stop("data with >3 dimensions are not allowed.")

    if(!all(is.na(group))) {
       if (!all(weights == 1) & opt$verbose > 0)
            cat("Warning: weights ignored in sm.density.compare\n")
       return(sm.density.compare(x, group, h, model, ...))
       }

    replace.na(opt, nbins, round((nobs > 500) * 8 * log(nobs) / ndim))
    rawdata <- list(nbins = opt$nbins, x = x, nobs = nobs, ndim = ndim)
    if(opt$nbins > 0) {
        if (!all(weights == 1) & opt$verbose>0)
            cat("Warning: weights overwritten by binning\n")
        bins    <- binning(x, nbins = opt$nbins)
        x       <- bins$x
        weights <- bins$x.freq
        nx      <- length(bins$x.freq)
        if(!all(is.na(opt$h.weights))) 
            stop("use of h.weights is incompatible with binning - set nbins = 0")
        }
    else
        nx <- nobs

    if(opt$positive) {
        replace.na(opt, delta, apply(as.matrix(x), 2, min))
        if ((ndim == 3 ) & (opt$verbose > 0)) 
           cat("the positive estimation is not available for 3 variables.\n")
        }   

    if(missing(h)){ 
        if(opt$positive) {
            xlog <- log(as.matrix(x) + outer(rep(1, nx), opt$delta))
            if (ndim == 1) xlog <- as.vector(xlog)
            h <- h.select(xlog, y = NA, weights = weights, ...)
            }
        else 
            h <- h.select(x = x, y = NA, weights = weights, ...)
    }
    
    if (opt$panel) 
       if (!require(rpanel)) {
       	  if (opt$verbose > 0) cat("The rpanel package is not available.\n")
          opt$panel <- FALSE
          }
    if (is.na(opt$band)) {
       if (model == "none") opt$band <- FALSE
          else              opt$band <- TRUE
       }
    if ((model == "none") && opt$band) opt$band <- FALSE

    if (ndim == 1) { 
        if(length(h)!=1) stop("length(h) != 1")
        replace.na(opt, xlab, x.name)
        replace.na(opt, ylab, "Probability density function")
        if (!opt$panel)
           est <- sm.density.1d(x, h, model, weights, rawdata, options = opt)
        else
           rp.density1(x, h, model, weights, rawdata, opt)
        }
    
    if (ndim == 2) {
        if(length(h) != 2) stop("length(h) != 2")
        dimn <- dimnames(x)[[2]]
        name.comp <- if (length(dimn) == 2) dimn
                     else outer(x.name, c("[1]","[2]"), paste, sep="")
        replace.na(opt, xlab, name.comp[1])
        replace.na(opt, ylab, name.comp[2])
        if (!opt$panel)
           est <- sm.density.2d(x, h, weights, rawdata, options = opt)
        else
           rp.density2(x, h, model, weights, rawdata, opt)
        }
    
    if (ndim == 3) {
        dimn <- dimnames(x)[[2]]
        name.comp <- if (length(dimn) == 3) dimn
                     else outer(x.name,c("[1]","[2]","[3]"),paste,sep="")
        replace.na(opt, xlab, name.comp[1])
        replace.na(opt, ylab, name.comp[2])
        replace.na(opt, zlab, name.comp[3])
        opt$nbins <- 0
        if (!opt$panel)
           est <- sm.density.3d(x, h = h, weights, rawdata, options = opt)
        else
           rp.density3(x, h, model, weights, rawdata, opt)
        }

    if (!opt$panel) {
       est$data <- list(x = x, nbins = opt$nbins, freq = weights)
       est$call <- match.call()
       invisible(est)
       }
    else
       invisible()
}

"sm.density.1d" <- function (x, h = hnorm(x, weights), model = "none", weights,
                                rawdata = list(x = x), options = list()) {

    absent <- function(x) missing(x) | any(is.na(x) | is.null(x))
    opt <- sm.options(options)
    replace.na(opt, display,    "line")
    replace.na(opt, col,        "black")
    replace.na(opt, col.band,   "cyan")
    replace.na(opt, col.points, "black")
    replace.na(opt, se,         FALSE)
    replace.na(opt, ngrid,      100)
    panel <- opt$panel
    band  <- opt$band
    hmult <- opt$hmult
    if (any(is.na(opt$h.weights)))
        replace.na(opt, h.weights, rep(1, length(x)))
    else band <- panel <- FALSE
    if (model == "none")
        band <- FALSE
    if (opt$add | opt$display %in% "none")
        panel <- FALSE
    a <- if (opt$positive) c(1 / opt$ngrid, max(x) * 1.05)
        else c(min(x) - diff(range(x)) / 4, max(x) + diff(range(x)) / 4)
    replace.na(opt, xlim, a)
    long.x <- rep(x, weights)
    a <- if (opt$positive)
            max(0.4/(quantile(long.x, 0.75) - quantile(long.x, 0.5)),
                0.4/(quantile(long.x, 0.5) - quantile(long.x, 0.25)))
        else 0.6/sqrt(wvar(x, weights))
    replace.na(opt, yht, a)
    replace.na(opt, ylim, c(0, opt$yht))
    replace.na(opt, ngrid, 100)
    if (!opt$add & !(opt$display %in% "none"))
        plot(opt$xlim, opt$ylim, type = "n", xlab = opt$xlab, ylab = opt$ylab)
    opt$band <- band
    opt$panel <- panel
    if (!(opt$display %in% "none"))
        est <- smplot.density(x, h, weights, rawdata, options = opt)
    if (all(!is.na(opt$eval.points))) {
        if (opt$positive)
            est <- sm.density.positive.1d(x, h, weights, options = opt)
        else est <- sm.density.eval.1d(x, h, weights, options = opt)
        }
    else if (opt$display %in% "none")
        est <- sm.density.eval.1d(x, h, weights = weights, options = opt)
    if (all(opt$h.weights == rep(1, length(x))) & opt$positive == FALSE) {
        se <- sqrt(dnorm(0, sd = sqrt(2))/(4 * sum(weights) *
            h))
        upper <- sqrt(est$estimate) + 2 * se
        lower <- pmax(sqrt(est$estimate) - 2 * se, 0)
        upper <- upper^2
        lower <- lower^2
        est$se <- rep(se, length(est$eval.points))
        est$upper <- upper
        est$lower <- lower
    }
    invisible(est)
}

"sm.density.2d" <- function (X, h = hnorm(X, weights), weights = rep(1, length(x)),
          rawdata = list(), options = list()) {
          	
    opt <- sm.options(options)

    x <- X[, 1]
    y <- X[, 2]
    replace.na(opt, display, "persp")
    if (opt$display == "contour") opt$display <- "slice"
    replace.na(opt, ngrid, 50)
    replace.na(opt, xlab, deparse(substitute(x)))
    replace.na(opt, ylab, deparse(substitute(y)))
    replace.na(opt, zlab, "Density function")
    if (any(is.na(opt$eval.points))) {
        replace.na(opt, xlim, range(X[, 1]))
        replace.na(opt, ylim, range(X[, 2]))
        replace.na(opt, eval.points,
            cbind(seq(opt$xlim[1], opt$xlim[2], length = opt$ngrid),
                  seq(opt$ylim[1], opt$ylim[2], length = opt$ngrid)))
        }
    else {
        replace.na(opt, xlim, range(opt$eval.points[, 1]))
        replace.na(opt, ylim, range(opt$eval.points[, 2]))
        }
    replace.na(opt, h.weights, rep(1, length(x)))
    hmult <- opt$hmult
    display <- opt$display
    if ((display == "rgl") & (!require(rgl))) {
        display <- "persp"
        cat("The rgl package is not available.\n")
        }
    if ((display == "rgl") & (!require(rpanel))) {
        display <- "persp"
        cat("The rpanel package is not available.\n")
        }
    surf.ids <- rep(NA, 2)

    if (!opt$eval.grid)
        est <- sm.density.eval.2d(x, y, h,
            xnew = opt$eval.points[,1], ynew = opt$eval.points[, 2],
            eval.type = "points", weights = weights, options = opt)
    else if (display == "none")
        est <- sm.density.eval.2d(x, y, h,
            xnew = opt$eval.points[,1], ynew = opt$eval.points[, 2],
            eval.type = "grid", weights = weights, options = opt)
    else if (display == "persp")
        est <- sm.persplot(x, y, h, weights, rawdata, options = opt)
    else if (display == "image")
        est <- sm.imageplot(x, y, h, weights, rawdata, options = opt)
    else if (display == "slice")
        est <- sm.sliceplot(x, y, h, weights, rawdata, options = opt)
    else if (display == "rgl")
        est <- sm.rglplot(x, y, h, weights, rawdata, options = opt)
    else
        stop("invalid setting for display.")

    if (all(opt$h.weights == rep(1, length(x)))) {
        se <- dnorm(0, sd = sqrt(2)) / sqrt(4 * sum(weights) * h[1] * h[2])
        upper <- sqrt(est$estimate) + 2 * se
        lower <- pmax(sqrt(est$estimate) - 2 * se, 0)
        upper <- upper^2
        lower <- lower^2
        est$se <- est$estimate - est$estimate + se
        est$upper <- upper
        est$lower <- lower
    }

    invisible(est)
}


"smplot.density" <- function (x, h, weights = rep(1, length(x)), rawdata = list(x = x),
    options = list()) {

    opt <- sm.options(options)
    if (opt$positive)
        est <- sm.density.positive.1d(x, h, weights = weights,
            options = opt)
    else {
        est <- sm.density.eval.1d(x, h, weights = weights, options = opt)
        if (opt$band)
            normdens.band(x, h, weights, options = opt)
        else if (!opt$add)
            polygon(c(par()$usr[1:2], par()$usr[2:1]), rep(c(par()$usr[3],
                par()$usr[4] * 0.999), c(2, 2)), col = 0, border = 0)
    }
    box()
    lines(est$eval.points, est$estimate, lty = opt$lty, col = opt$col, lwd = opt$lwd)
    if (opt$rugplot && !opt$add)
        rug(jitter(rawdata$x, amount = 0), 0.015)
    if ((opt$se | opt$display %in% "se") & (!opt$band) & 
                       all(opt$h.weights == rep(1, length(x)))) {
        se <- sqrt(dnorm(0, sd = sqrt(2))/(4 * h * sum(weights)))
        upper <- sqrt(est$estimate) + 2 * se
        lower <- pmax(sqrt(est$estimate) - 2 * se, 0)
        upper <- upper^2
        lower <- lower^2
        lines(est$eval.points, upper, lty = 3, col = opt$col)
        lines(est$eval.points, lower, lty = 3, col = opt$col)
    }
    invisible(est)
}

"normdens.band" <- function (x, h, weights = rep(1, length(x)), options = list()) {

    opt <- sm.options(options)
    xlim <- opt$xlim
    yht <- opt$yht
    ngrid <- opt$ngrid
    x.points <- seq(xlim[1], xlim[2], length = ngrid)
    xbar <- wmean(x, weights)
    sx <- sqrt(wvar(x, weights))
    hm <- h * opt$hmult
    dmean <- dnorm(x.points, xbar, sqrt(sx^2 + hm^2))
    dvar <- (dnorm(0, 0, sqrt(2 * hm^2)) * dnorm(x.points, xbar,
        sqrt(sx^2 + 0.5 * hm^2)) - (dmean)^2)/sum(weights)
    upper <- pmin(dmean + 2 * sqrt(dvar), par()$usr[4])
    lower <- pmax(0, dmean - 2 * sqrt(dvar))
    polygon(c(par()$usr[1:2], par()$usr[2:1]), rep(c(par()$usr[3],
        par()$usr[4] * 0.999), c(2, 2)), col = 0, border = FALSE)
    polygon(c(x.points, rev(x.points)), c(upper, rev(lower)),
            col = opt$col.band, border = FALSE)
}

"sm.density.compare" <- function (x, group, h, model = "none", ...) {

    if (!is.vector(x))
       stop("sm.density.compare can handle only 1-d data") 
    opt <- sm.options(list(...))
    replace.na(opt, ngrid, 50)
    replace.na(opt, display, "line")
    replace.na(opt, xlab, deparse(substitute(x)))
    replace.na(opt, ylab, "Density")
    replace.na(opt, xlim, c(min(x) - diff(range(x))/4, max(x) +
                            diff(range(x))/4))
    replace.na(opt, eval.points,
        seq(opt$xlim[1], opt$xlim[2], length=opt$ngrid))
    if (is.na(opt$band)) {
       if (model == "none") opt$band <- FALSE
          else              opt$band <- TRUE
       }
    if ((model == "none") && opt$band) opt$band <- FALSE
    band <- opt$band
    ngrid <- opt$ngrid
    xlim <- opt$xlim
    nboot <- opt$nboot
    y <- x
    if (is.na(opt$test)) {
       if (model == "none") opt$test <- FALSE
          else              opt$test <- TRUE
       }
    if ((model == "none") && opt$test) opt$test <- FALSE
    test <- opt$test
    if (opt$display %in% "none")
        band <- FALSE
    fact <- factor(group)
    fact.levels <- levels(fact)
    nlev <- length(fact.levels)
    ni <- table(fact)
    if (band & (nlev > 2)) {
        cat("Reference band available to compare two groups only.", "\n")
        band <- FALSE
    }
    if (length(opt$lty) < nlev) opt$lty <- 1:nlev
    if (length(opt$col) < nlev) opt$col <- 2:(nlev + 1)
    if (missing(h))
        h <- h.select(x, y = NA, group = group, ...)
    opt$band <- band
    opt$test <- test
    estimate <- matrix(0, ncol = opt$ngrid, nrow = nlev)
    se <- matrix(0, ncol = opt$ngrid, nrow = nlev)
    for (i in 1:nlev) {
        sm <- sm.density(y[fact == fact.levels[i]], h = h, display = "none",
            eval.points = opt$eval.points)
        estimate[i, ] <- sm$estimate
        se[i, ] <- sm$se
    }
    eval.points <- sm$eval.points
    if (!(opt$display %in% "none" | band)) {
        replace.na(opt, yht, 1.1 * max(as.vector(estimate)))
        replace.na(opt, ylim, c(0, opt$yht))
        plot(xlim, opt$ylim, xlab = opt$xlab, ylab = opt$ylab, type = "n")
        for (i in 1:nlev) lines(eval.points, estimate[i, ],
	          lty = opt$lty[i], col = opt$col[i], lwd = opt$lwd)
    }
    est <- NULL
    p <- NULL
    if (model == "equal" & test) {
        if (nlev == 2) {
            ts <- sum((estimate[1, ] - estimate[2, ])^2)
        }
        else {
            sm.mean <- sm.density(y, h = h, xlim = opt$xlim,
                ngrid = opt$ngrid, display = "none")$estimate
            ts <- 0
            for (i in 1:nlev) ts <- ts + ni[i] *
                sum((estimate[i,] - sm.mean)^2)
        }
        p <- 0
        est.star <- matrix(0, ncol = opt$ngrid, nrow = nlev)
        for (iboot in 1:nboot) {
            ind <- (1:length(y))
            for (i in 1:nlev) {
                indi <- sample((1:length(ind)), ni[i])
                est.star[i, ] <- sm.density(y[ind[indi]], h = h,
                  ngrid = opt$ngrid, xlim = opt$xlim, display = "none")$estimate
                ind <- ind[-indi]
            }
            if (nlev == 2) {
                ts.star <- sum((est.star[1, ] - est.star[2, ])^2)
            }
            else {
                sm.mean <- sm.density(y, h = h, xlim = opt$xlim,
                  ngrid = opt$ngrid, display = "none")$estimate
                ts.star <- 0
                for (i in 1:nlev) {
                  ts.star <- ts.star + ni[i] * sum((est.star[i,] - sm.mean)^2)
                  }
                }
            if (ts.star > ts)
                p <- p + 1
            if (opt$verbose > 1) {
                cat(iboot)
                cat(" ")
            }
        }
        p <- p/nboot
        cat("\nTest of equal densities:  p-value = ", round(p,3), "\n")
        est <- list(p = p, estimaate = estimate, eval.points = eval.points, h = h)
    }
    if (model == "equal" & band) {
        av <- (sqrt(estimate[1, ]) + sqrt(estimate[2, ]))/2
        se <- sqrt(se[1, ]^2 + se[2, ]^2)
        upper <- (av + se)^2
        lower <- pmax(av - se, 0)^2
        replace.na(opt, yht, 1.1 * max(as.vector(estimate), upper))
        replace.na(opt, ylim, c(0, opt$yht))
        plot(xlim, opt$ylim, xlab = opt$xlab, ylab = opt$ylab, type = "n")
        polygon(c(eval.points, rev(eval.points)), c(upper, rev(lower)),
            col = "cyan", border = 0)
        lines(eval.points, estimate[1, ], lty = opt$lty[1], col = opt$col[1], lwd = opt$lwd)
        lines(eval.points, estimate[2, ], lty = opt$lty[2], col = opt$col[2], lwd = opt$lwd)
        est <- list(p = p, estimate = estimate, eval.points = eval.points, 
                    upper = upper, lower = lower, h = h)
    }
    invisible(est)
}

"sm.density.eval.1d" <- function (x, h, weights = rep(1, n), options = list()) {
	
    opt <- sm.options(options)
    replace.na(opt, h.weights, rep(1, length(x)))
    replace.na(opt, xlim, c(min(x) - diff(range(x))/4, max(x) +
        diff(range(x))/4))
    replace.na(opt, ngrid, 100)
    hmult <- opt$hmult
    h.weights <- opt$h.weights
    xlim <- opt$xlim
    ngrid <- opt$ngrid
    replace.na(opt, eval.points, seq(xlim[1], xlim[2], length = ngrid))
    xnew <- opt$eval.points
    n <- length(x)
    neval <- length(xnew)
    W <- matrix(rep(xnew, rep(n, neval)), ncol = n, byrow = TRUE)
    W <- W - matrix(rep(x, neval), ncol = n, byrow = TRUE)
    W1 <- matrix(rep(h.weights, neval), ncol = n, byrow = TRUE)
    W <- exp(-0.5 * (W/(hmult * h * W1))^2)/W1
    est <- W %*% weights/(sum(weights) * sqrt(2 * pi) * hmult * h)
    invisible(list(eval.points = xnew, estimate = as.vector(est),
        h = h * hmult, h.weights = h.weights, weights = weights))
}

"sm.density.eval.2d" <- function (x, y, h, xnew, ynew, eval.type = "points", weights = rep(1, n),
          options = list()) {

    opt <- sm.options(options)
    replace.na(opt, xlim, range(x))
    replace.na(opt, ylim, range(y))
    replace.na(opt, ngrid, 50)
    replace.na(opt, h.weights, rep(1, length(x)))
    if (missing(xnew))
        xnew <- seq(opt$xlim[1], opt$xlim[2], length = opt$ngrid)
    if (missing(ynew))
        ynew <- seq(opt$ylim[1], opt$ylim[2], length = opt$ngrid)
    n <- length(x)
    nnew <- length(xnew)
    h.weights <- opt$h.weights
    hmult <- opt$hmult
    W1 <- matrix(rep(xnew, rep(n, nnew)), ncol = n, byrow = TRUE)
    W1 <- W1 - matrix(rep(x, nnew), ncol = n, byrow = TRUE)
    W2 <- matrix(rep(h.weights, nnew), ncol = n, byrow = TRUE)
    Wx <- exp(-0.5 * (W1/(hmult * h[1] * W2))^2)/W2
    W1 <- matrix(rep(ynew, rep(n, nnew)), ncol = n, byrow = TRUE)
    W1 <- W1 - matrix(rep(y, nnew), ncol = n, byrow = TRUE)
    Wy <- exp(-0.5 * (W1/(opt$hmult * h[2] * W2))^2)/W2
    if (eval.type == "points")
        est <- as.vector(((Wx * Wy) %*% weights)/(sum(weights) *
            2 * pi * h[1] * h[2] * hmult^2))
    else est <- (Wx %*% (weights * t(Wy)))/(sum(weights) * 2 *
        pi * h[1] * h[2] * hmult^2)
    invisible(list(eval.points = cbind(xnew, ynew), estimate = est,
        h = h * hmult, h.weights = h.weights, weights = weights))
}

"sm.density.positive.1d" <- function (x, h, weights, options = list()) {

    opt <- sm.options(options)
    if (min(x) <= 0 & opt$verbose>0)
        cat("Warning: some data are not positive\n")
    delta <- opt$delta
    replace.na(opt, ngrid, 100)
    replace.na(opt, xlim, c(min(x), max(x)))
    if (min(opt$xlim) < 0 & opt$verbose>0)
        cat("Warning: xlim<0 with positive=TRUE \n")
    if (missing(h))
        h <- hnorm(log(x + delta), weights)
    ngrid <- opt$ngrid
    ev.pt <- opt$eval.points
    if (any(is.na(ev.pt))) {
        a     <- log(opt$xlim + 1/ngrid)
        ev.pt <- exp(seq(min(a), max(a), length=opt$ngrid))
    }
    opt$eval.points <- log(ev.pt + delta)
    f <- sm.density.eval.1d(log(x + delta), h = h, weights = weights,
        options = opt)
    est <- f$estimate/(ev.pt + delta)
    est[is.na(est)] <- 0
    list(eval.points = ev.pt, estimate = as.vector(est), h = h)
}

"sm.density.positive.2d" <- function (X, h = c(hnorm(log(X[, 1] + delta[1]), weights), 
            hnorm(log(X[,2] + delta[2]), weights)), eval.type = "points", 
            weights = rep(1, nrow(X)), options = list()) {

    opt <- sm.options(options)
    replace.na(opt, ngrid, 50)
    replace.na(opt, delta, apply(X, 2, min))
    if (min(X) <= 0 & opt$verbose > 0)
        cat("Warning: some data are not positive\n")
    if (dim(X)[2] != 2)
        cat("parameter X must be a two-column matrix\n")
    x1 <- X[, 1]
    x2 <- X[, 2]
    delta <- opt$delta
    replace.na(opt, xlim, range(x1))
    replace.na(opt, ylim, range(x2))
    replace.na(opt, ngrid, 50)
    xlim <- opt$xlim
    ylim <- opt$ylim
    ngrid <- opt$ngrid
    ax <- log(xlim + 1/ngrid)
    ay <- log(ylim + 1/ngrid)
    eval1 <- exp(seq(ax[1], ax[2], length = ngrid)) - 1/ngrid
    eval2 <- exp(seq(ay[1], ay[2], length = ngrid)) - 1/ngrid
    replace.na(opt, eval.points, cbind(eval1, eval2))
    eval1 <- opt$eval.points[, 1]
    eval2 <- opt$eval.points[, 2]
    pdf <- sm.density.eval.2d(log(x1 + delta[1]), log(x2 + delta[2]),
        h = h, xnew = log(eval1 + delta[1]), ynew = log(eval2 +
            delta[2]), eval.type = eval.type, weights = weights)
    if (eval.type == "points")
        est <- pdf$estimate/((eval1 + delta[1]) * (eval2 + delta[2]))
    else est <- pdf$estimate/outer(eval1 + delta[1], eval2 +
        delta[2])
    invisible(list(x1 = eval1, x2 = eval2, estimate = est, h = h))
}

"sm.density.positive.grid" <- function (X, 
    h = c(hnorm(log(X[, 1] + delta[1])), hnorm(log(X[, 2] + delta[2]))), 
    weights=NA, options=list()) {

    f <- sm.density.positive.2d(X, h, eval.type = "grid", 
                          weights=weights, options=options)
    invisible(list(eval.points = cbind(f$x1, f$x2), estimate = f$est,
        h = h))
}

"sm.imageplot" <- function (x, y, h, weights, rawdata, options = list()) {

    opt <- sm.options(options)

    ngrid <- opt$ngrid
    xlim <- opt$xlim
    ylim <- opt$ylim
    xgrid <- opt$eval.points[,1]
    ygrid <- opt$eval.points[,2]
    if (!opt$positive)
        dgrid <- sm.density.eval.2d(x, y, h, xgrid, ygrid,
	    eval.type = "grid", weights, opt)$estimate
    else {
        f <- sm.density.positive.grid(cbind(x, y), h, weights=weights, 
                options=opt)
        xgrid <- f$eval.points[, 1]
        ygrid <- f$eval.points[, 2]
        dgrid <- f$estimate
        }
    image(xgrid, ygrid, dgrid,
          xlab = opt$xlab, ylab = opt$ylab, xlim = xlim, ylim = ylim,
          add = opt$add, col = opt$col.palette)
    invisible(list(eval.points = cbind(xgrid, ygrid), estimate = dgrid,
        h = h * opt$hmult, h.weights = opt$h.weights, weights = weights))
    }


"sm.persplot" <- function (x, y, h = hnorm(cbind(x, y), weights), weights, rawdata = list(),
    options = list()) {

    opt <- sm.options(options)

    ngrid <- opt$ngrid
    xlim <- opt$xlim
    ylim <- opt$ylim
    xgrid <- opt$eval.points[,1]
    ygrid <- opt$eval.points[,2]
    if (!opt$positive)
        dgrid <- sm.density.eval.2d(x, y, h, xgrid, ygrid, eval.type = "grid",
            weights, options = opt)$estimate
    else {
        f <- sm.density.positive.grid(cbind(x, y), h, weights=weights, 
                       options=opt)
        xgrid <- f$eval.points[, 1]
        ygrid <- f$eval.points[, 2]
        dgrid <- f$estimate
    }
    zlim <- replace.na(opt, zlim, c(0, max(dgrid, na.rm = TRUE)))
    if (is.na(opt$col)) opt$col <- "green"
    persp(xgrid, ygrid, dgrid,
        xlab = opt$xlab, ylab = opt$ylab, zlab = opt$zlab,
        xlim = xlim, ylim = ylim, zlim = opt$zlim,
        theta = opt$theta, phi = opt$phi,
        ticktype = "detailed", col = opt$col, d = 4)
    invisible(list(eval.points = cbind(xgrid, ygrid), estimate = dgrid,
        h = h * opt$hmult, h.weights = opt$h.weights, weights = weights))
}


"sm.sliceplot" <- function (x, y, h, weights, rawdata = list(), options = list()) {

    opt <- sm.options(options)

    ngrid <- opt$ngrid
    xlim  <- opt$xlim
    ylim  <- opt$ylim
    xgrid <- opt$eval.points[,1]
    ygrid <- opt$eval.points[,2]

    if (!opt$add) {
        plot(x, y, xlim = opt$xlim, ylim = opt$ylim, xlab = opt$xlab,
            ylab = opt$ylab, type = "n")
        points(rawdata$x[, 1], rawdata$x[, 2], col = opt$col,
            pch = opt$pch, cex = 2/log(rawdata$nobs))
    }

    if (opt$positive)
        f <- sm.density.positive.grid(cbind(x, y), h, weights=weights, 
                      options=opt)
    else f <- sm.density.eval.2d(x, y, h, xgrid, ygrid,
        eval.type = "grid", weights = weights, options = opt)
    dgrid <- f$estimate
    xgrid <- f$eval.points[, 1]
    ygrid <- f$eval.points[, 2]
    if (opt$positive) {
        opt$eval.points <- cbind(x, y)
        dobs <- sm.density.positive.2d(cbind(x, y), h, eval.type = "points",
            weights = weights, options = opt)$estimate
    }
    else dobs <- sm.density.eval.2d(x, y, h, xnew = x, ynew = y,
        weights = weights)$estimate
    props <- opt$props
    hts <- quantile(rep(dobs, weights), prob = (100 - props)/100)
    for (i in 1:length(props)) {
        scale <- props[i]/hts[i]
        contour(xgrid, ygrid, dgrid * scale, level = hts[i] *
            scale, add = TRUE, lty = opt$lty, col = opt$col)
    }
    invisible(list(eval.points = cbind(xgrid, ygrid), estimate = dgrid,
        h = h * opt$hmult, h.weights = opt$h.weights, weights = weights))
}

"sm.rglplot" <- function (x, y, h, weights, rawdata, options = list()) {

    opt <- sm.options(options)

    ngrid <- opt$ngrid
    xlim  <- opt$xlim
    ylim  <- opt$ylim
    xgrid <- opt$eval.points[, 1]
    ygrid <- opt$eval.points[, 2]
    if (!opt$positive)
        dgrid <- sm.density.eval.2d(x, y, h, xgrid, ygrid,
	                          eval.type = "grid", weights, opt)$estimate
    else {
        f     <- sm.density.positive.grid(cbind(x, y), h, weights = weights, options = opt)
        xgrid <- f$eval.points[, 1]
        ygrid <- f$eval.points[, 2]
        dgrid <- f$estimate
        }
     
     if (!opt$add) {
     	if (any(is.na(opt$zlim))) opt$zlim <- c(0, max(dgrid * 1.5))
        opt$scaling <- rp.plot3d(xgrid, dgrid, ygrid, type = "n",
                         xlab = opt$xlab, ylab = opt$zlab, zlab = opt$ylab,
                         xlim = opt$xlim, ylim = opt$zlim, zlim = opt$ylim,
                         size = opt$size, col = opt$col.points)
        }
     surf.ids <- sm.surface3d(cbind(xgrid, ygrid), dgrid, opt$scaling, 
                   col = opt$col, col.mesh = opt$col.mesh, alpha = opt$alpha, 
                   lit = opt$lit)
    invisible(list(eval.points = cbind(xgrid, ygrid), estimate = dgrid,
        h = h * opt$hmult, h.weights = opt$h.weights, weights = weights, 
        scaling = opt$scaling, surf.ids = surf.ids))
    }

"sm.density.3d" <- function(x, h = hnorm(x, weights), 
                             weights = rep(1, length(x)), rawdata = list(), 
                             options = list()) {
          	
    opt <- sm.options(options)

    replace.na(opt, ngrid, 20)
    replace.na(opt, xlab, deparse(substitute(x)))
    replace.na(opt, ylab, deparse(substitute(y)))
    replace.na(opt, zlab, deparse(substitute(z)))
    replace.na(opt, display, "rgl")
    if (any(is.na(opt$col)) | length(opt$col) != length(opt$props))
       opt$col <- topo.colors(length(opt$props))
    if (length(opt$alpha) != length(opt$props))
       opt$alpha <- seq(1, 0.5, length = length(opt$props))
    if (any(is.na(opt$eval.points))) {
        replace.na(opt, xlim, range(x[, 1]))
        replace.na(opt, ylim, range(x[, 2]))
        replace.na(opt, zlim, range(x[, 3]))
        evp <- cbind(seq(opt$xlim[1], opt$xlim[2], length = opt$ngrid),
                     seq(opt$ylim[1], opt$ylim[2], length = opt$ngrid),
                     seq(opt$zlim[1], opt$zlim[2], length = opt$ngrid))
        replace.na(opt, eval.points, evp)
        }
    else {
        replace.na(opt, xlim, range(opt$eval.points[, 1]))
        replace.na(opt, ylim, range(opt$eval.points[, 2]))
        replace.na(opt, zlim, range(opt$eval.points[, 3]))
        }
    replace.na(opt, h.weights, rep(1, nrow(x)))
    hmult <- opt$hmult
    display <- opt$display

   est <- sm.density.eval.3d(x, h, opt$eval.points, eval.type = "grid", 
                  weights = weights, rawdata = rawdata, options = opt)

   invisible(est)
   }

"sm.density.eval.3d" <- function (x, h, eval.points, eval.type = "points", 
                         weights = rep(1, nrow(x)), rawdata = list(), options = list()) {

    opt <- sm.options(options)
    replace.na(opt, xlim, range(x[, 1]))
    replace.na(opt, ylim, range(x[, 2]))
    replace.na(opt, zlim, range(x[, 3]))
    replace.na(opt, ngrid, 20)
    replace.na(opt, h.weights, rep(1, nrow(x)))
    n         <- nrow(x)
    nnew      <- nrow(eval.points)
    h.weights <- opt$h.weights
    hmult     <- opt$hmult
    result    <- list(eval.points = eval.points,
                      h = h * hmult, h.weights = h.weights, weights = weights)
    surf.ids <- NA
    
    Wh <- matrix(rep(h.weights, nnew), ncol = n, byrow = TRUE)
    W1 <- matrix(rep(eval.points[, 1], rep(n, nnew)), ncol = n, byrow = TRUE)
    W1 <- W1 - matrix(rep(x[, 1], nnew), ncol = n, byrow = TRUE)
    W1 <- exp(-0.5 * (W1/(hmult * h[1] * Wh))^2)/Wh
    W2 <- matrix(rep(eval.points[, 2], rep(n, nnew)), ncol = n, byrow = TRUE)
    W2 <- W2 - matrix(rep(x[, 2], nnew), ncol = n, byrow = TRUE)
    W2 <- exp(-0.5 * (W2/(hmult * h[2] * Wh))^2)/Wh
    W3 <- matrix(rep(eval.points[, 3], rep(n, nnew)), ncol = n, byrow = TRUE)
    W3 <- W3 - matrix(rep(x[, 3], nnew), ncol = n, byrow = TRUE)
    W3 <- exp(-0.5 * (W3/(hmult * h[3] * Wh))^2)/Wh
    
    if (eval.type == "points")
       est <- as.vector(((W1 * W2 * W3) %*% weights) / (sum(weights) *
            (2 * pi)^1.5 * h[1] * h[2] * h[3] * hmult^3))
    else {
       est <- sm.density.eval.3d(x, h, x, eval.type = "points", 
                       weights = weights, options = opt)$estimate
       levels <- quantile(est, opt$props / 100)
       est <- apply(W3, 1, function(x) (W1 %*% (weights * x * t(W2))) /
                      (sum(weights) * (2 * pi)^1.5 * h[1] * h[2] * h[3] * hmult^3))
       est <- array(c(est), dim = rep(opt$ngrid, 3))
       if (opt$display != "none") {
          if (require(rpanel) & require(rgl) & require(misc3d)) {
             struct <- contour3d(est, levels, 
                           eval.points[, 1], eval.points[, 2], eval.points[, 3], 
                           engine = "none")
             if (!opt$add) {
             	opt$scaling <- rp.plot3d(rawdata$x[, 1], rawdata$x[, 2], rawdata$x[, 3],
                                  xlab = opt$xlab, ylab = opt$ylab, zlab = opt$zlab,
                                  col = opt$col.points, size = opt$size)
                result$scaling <- opt$scaling
                }
             surf.ids <- integer(0)
             for (i in 1:length(opt$props)) {
             	if (length(opt$props) > 1) strct <- struct[[i]] else strct <- struct
                trngs.x <- c(t(cbind(strct$v1[, 1], strct$v2[, 1],strct$v3[, 1])))
                trngs.y <- c(t(cbind(strct$v1[, 2], strct$v2[, 2],strct$v3[, 2])))
                trngs.z <- c(t(cbind(strct$v1[, 3], strct$v2[, 3],strct$v3[, 3])))
                a <- opt$scaling(trngs.x, trngs.y, trngs.z)
                surf.ids <- c(surf.ids, 
                      triangles3d(a$x, a$y, a$z, col = opt$col[i], alpha = opt$alpha[i]))
                }
             }
          else if (opt$verbose > 0) cat("at least one of the rpanel, rgl or misc3d packages",
                   " is not available.\n")
          }
       }
    result$estimate <- est
    result$surf.ids <- surf.ids
    invisible(result)
    }


"nise" <- function (y, ...) {
    n <- length(y)
    opt <- sm.options(list(...))
    replace.na(opt, nbins, round((n > 500) * 8 * log(n)))
    replace.na(opt, hmult, 1)
    if (opt$nbins > 0) {
        bins <- binning(y, nbins = opt$nbins)
        y <- bins$x
        weights <- bins$x.freq
    }
    else weights <- rep(1, n)
    y <- (y - wmean(y, weights)) / sqrt(wvar(y, weights))
    h <- hnorm(y) * opt$hmult
    result <- dnorm(0, sd = sqrt(2 + 2 * h^2))
    result <- result - 2 * sm.density(y, h = sqrt(1 + 2 * h^2),
        eval.points = 0, display = "none", weights = weights,
        nbins = 0)$estimate
    result <- result + wmean(sm.density(y, h = sqrt(2) * h, eval.points = y,
        display = "none", weights = weights, nbins = 0)$estimate,
        weights)
    result
    }

"nmise" <- function (sd, n, h) {
    dnorm(0, sd = sqrt(2) * h)/n + (1 - 1/n) * dnorm(0, sd = sqrt(2 *
        (sd^2 + h^2))) - 2 * dnorm(0, sd = sqrt(2 * sd^2 + h^2)) +
        dnorm(0, sd = sqrt(2) * sd)
    }

"nnbr" <- function (x, k) {
    if (isMatrix(x)) {
        ndim <- 2
        n <- nrow(x)
        }
    else {
        ndim <- 1
        n <- length(x)
        }
    knn <- vector("numeric", n)
    if (ndim == 1) {
        for (i in 1:length(x)) knn[i] <- sort(abs(x - x[i]))[k +
            1]}
    if (ndim == 2) {
        for (i in 1:length(x[, 1])) knn[i] <- sort(sqrt(((x[,
            1] - x[i, 1])^2)/var(x[, 1]) + ((x[, 2] - x[i, 2])^2)/var(x[,
            2])))[k + 1]
        }
    knn
    }
