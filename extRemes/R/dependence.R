atdf <- function(x, u, lag.max=NULL, type=c("all", "rho", "rhobar"), plot=TRUE, na.action=na.fail, ...) {
   out <- list()
   out$call <- match.call()
   if(!is.numeric(x)) stop("atdf: x must be numeric.")
   if(u <= 0 || u >= 1) stop("atdf: u must be between 0 and 1 (non-inclusive).")
   type <- match.arg(type)
   out$type <- type
   out$series <- deparse(substitute(x))
   x <- na.action(as.ts(x))
   n <- length(x)
   if(is.null(lag.max)) lag.max <- floor(10*(log10(n)))
   lag.max <- min(lag.max, n - 1)
   if(lag.max < 0) stop("atdf: lag.max must be at least 0.")
   lag <- 0:lag.max
   out$lag <- lag
   if(type=="all") {
	res <- matrix(NA, lag.max+1, 2)
	colnames(res) <- c("rho", "rhobar")
   } else res <- numeric(lag.max)+NA
   xun <- sort(x)[floor(n*u)]
   for(h in lag) {
	hold <- cbind(x[(1+h):n], x[1:(n - h)])
	hold <- apply(hold,1,min) > xun
	if(type=="all" || type=="rho") {
	   if(type=="all") res[h+1,1] <- sum(hold)/(n*(1-u))
	   else res[h+1] <- sum(hold)/(n*(1-u))
	}
	if(type=="all" || type=="rhobar") {
	   if(type=="all") res[h+1,2] <- -1 + 2*log(1-u)/log(sum(hold)/n)
	   else res[h+1] <- -1 + 2*log(1-u)/log(sum(hold)/n)
	}
   } # end of for 'h' loop.
   out$atdf <- res
   class(out) <- "atdf"
   if(plot) {
	plot(out, ...)
	invisible(out)
   } else return(out)
} # end of 'atdf' function.

plot.atdf <- function(x, type=NULL, ...) {
   if(is.null(type)) type <- x$type
   if(type=="all") {
	if(x$type != "all") stop("plot.atdf: invalid type argument.")
	x1 <- x$atdf[,1]
	x2 <- x$atdf[,2]
	par(mfrow=c(2,1), mar=c(4.1, 4.1, 2.1, 1.1))
   } else if(x$type=="all") {
	if(type=="rho") x1 <- x2 <- x$atdf[,1]
	else x1 <- x2 <- x$atdf[,2]
   } else x1 <- x2 <- x$atdf
   args <- list(...)
   if(is.null(args$main)) {
      if(is.null(args$ylab) && is.null(args$xlab)) {
         if(type=="all" ) plot(x$lag, x1, type="h", ylim=c(0,1), main=paste("auto-tail dependence function\n", x$series, sep=""), ylab="rho", xlab="", ...)
	 else if(type=="chi") plot(x$lag, x1, type="h", ylim=c(0,1), main=paste("auto-tail dependence function", x$series, sep=""), ylab="rho", xlab="lag", ...)
         if(type=="all") plot(x$lag, x2, type="h", ylim=c(-1,1), ylab="rhobar", xlab="lag", ...)
         else if(type=="rhobar") plot(x$lag, x2, type="h", ylim=c(-1,1), main=paste("auto-tail dependence function", x$series, sep=""), ylab="rhobar", xlab="lag", ...)
      } else if(is.null(args$ylab)) {
         if(type=="all" || type=="rho") plot(x$lag, x1, type="h", ylim=c(0,1), main=paste("auto-tail dependence function", x$series, sep=""), ylab="rho", ...)
         if(type=="all") plot(x$lag, x2, type="h", ylim=c(-1,1), ylab="rhobar", ...)
	 else if(type=="rhobar") plot(x$lag, x2, type="h", ylim=c(-1,1), main=paste("auto-tail dependence function", x$series, sep=""), ylab="rhobar", ...)
      } else if(is.null(args$xlab)) {
         if(type=="all") plot(x$lag, x1, type="h", main=paste("auto-tail dependence function", x$series, sep=""), ylim=c(0,1), xlab="", ...)
	 else if(type=="rho") plot(x$lag, x1, type="h", main=paste("auto-tail dependence function", x$series, sep=""), ylim=c(0,1), xlab="lag", ...)
         if(type=="all") plot(x$lag, x2, type="h", ylim=c(-1,1), xlab="lag", ...)
	 else if(type=="rhobar") plot(x$lag, x2, type="h", ylim=c(-1,1), main=paste("auto-tail dependence function", x$series, sep=""), xlab="lag", ...)
      } else {
         if(type=="all" || type=="rho") plot(x$lag, x1, type="h", ylim=c(0,1), main=paste("auto-tail dependence function", x$series, sep=""), ...)
         if(type=="all" || type=="rhobar") plot(x$lag, x2, type="h", ylim=c(-1,1), main=paste("auto-tail dependence function", x$series, sep=""), ...)
      }
   } else {
      if(is.null(args$ylab) && is.null(args$xlab)) {
         if(type=="all") plot(x$lag, x1, type="h", ylim=c(0,1), ylab="rho", xlab="", ...)
	 else if(type=="rho") plot(x$lag, x1, type="h", ylim=c(0,1), ylab="rho", xlab="lag", ...)
         if(type=="all") plot(x$lag, x2, type="h", ylim=c(-1,1), ylab="rhobar", xlab="lag", ...)
         else if(type=="rhobar") plot(x$lag, x2, type="h", ylim=c(-1,1), ylab="rhobar", xlab="lag", ...)
      } else if(is.null(args$ylab)) {
         if(type=="all" || type=="rho") plot(x$lag, x1, type="h", ylim=c(0,1), ylab="rho", ...)
         if(type=="all" || type=="rhobar") plot(x$lag, x2, type="h", ylim=c(-1,1), ylab="rhobar", ...)
      } else if(is.null(args$xlab)) {
         if(type=="all") plot(x$lag, x1, type="h", ylim=c(0,1), xlab="", ...)
	 else if(type=="rho") plot(x$lag, x1, type="h", ylim=c(0,1), xlab="lag", ...)
         if(type=="all" || type=="rhobar") plot(x$lag, x2, type="h", ylim=c(-1,1), xlab="lag", ...)
      } else {
         if(type=="all" || type=="rho") plot(x$lag, x1, type="h", ylim=c(0,1), ...)
         if(type=="all" || type=="rhobar") plot(x$lag, x2, type="h", ylim=c(-1,1), ...)
      }
   } # end of if else 'main' stmts.
   invisible()
} # end of 'plot.atdf' function.

taildep <- function(x, y, u, type=c("all", "chi", "chibar"), na.rm=FALSE) {
    type <- match.arg(type)
    if(missing(y)) {
	if((is.matrix(x) || is.data.frame(x) || is.matrix(x)) && dim(x)[2] == 2) y <- x[,2]
	else stop("taildep: must supply y or else x must be a two column matrix or data frame object.")
	x <- x[,1]
    }
    n <- length(x)
    if(n != length(y)) stop("taildep: y must have same length as x.")
    if(na.rm) {
	good <- !is.na(x) && !is.na(y)
	x <- x[good]
	y <- y[good]
	n <- length(x)
    }
    xun <- sort(x)[floor(n*u)]
    yun <- sort(y)[floor(n*u)]
    id <- (x > xun) & (y > yun)
    if(type=="all" || type=="chi") chi <- sum(id, na.rm=TRUE)/(n*(1-u))
    if(type=="all" || type=="chibar") chibar <- 2*log(1 - u)/log(mean(id)) - 1
    if(type=="all") {
	res <- c(chi, chibar)
	names(res) <- c("chi", "chibar")
    } else if(type=="chi") {
	res <- chi
	names(res) <- "chi"
    } else {
	res <- chibar
	names(res) <- "chibar"
    }
    return(res)
} # end of 'taildep' function.

taildep.test <- function(x, y, cthresh=-0.5, trans="relative.rank", na.action=na.fail, ...) {
    out <- list()
    out$call <- match.call()
    out$data.name <- c(deparse(substitute(x)),deparse(substitute(y)))
    out$method <- "Reiss-Thomas (13.35)"
    out$transformation <- trans
    PARAMETER <- c(cthresh, NA, NA, NA)
    names(PARAMETER) <- c("threshold", "m", "n", "exceedance rate (%)")
    if(missing(y)) {
        if((is.matrix(x) || is.data.frame(x)) && dim(x)[2] == 2) y <- x[,2]
        else stop("taildep.test: must supply y or else x must be a two column matrix or data frame object.")
        x <- x[,1]
    }
    tmp <- cbind(x,y)
    tmp <- na.action(tmp)
    n <- dim(tmp)[1]
    x <- tmp[,1]
    y <- tmp[,2]
    out$trans.args <- list(...)
    u <- do.call(trans, c(list(x=x), list(...)))
    v <- do.call(trans, c(list(x=y), list(...)))
    u <- u - 1
    v <- v - 1
    ci <- (u + v)
    out$c.full <- ci
    ci <- ci[ci > cthresh]
    m <- length(ci)
    PARAMETER[2:4] <- c(m, n, round(m/n, digits=3)*100)
    ci <- log(ci/cthresh)
    STATISTIC <- -(sum(ci,na.rm=TRUE) + m)/sqrt(m)
    names(STATISTIC) <- "statistic"
    out$statistic <- STATISTIC
    out$alternative <- "less"
    PVAL <- pnorm(STATISTIC, lower.tail=TRUE)
    out$parameter <- PARAMETER
    out$p.value <- PVAL
    class(out) <- "htest"
    return(out)
} # end of 'taildep.test' function.

relative.rank <- function(x, div="n", ...) {
    if(div=="n") n <- length(x)
    else if(is.element(div, c("n+1","n + 1", "n+ 1", "n +1"))) n <- length(x)+1
    else stop("relative.rank: div must be one of n or n+1.")
    x <- (rank(x, ...)/n)
    return(x)
} # end of 'relative.rank' funciton.

extremalindex <- function(x, threshold, method=c("intervals","runs"), run.length=1, na.action=na.fail, ...) {

    method <- tolower(method)
    method <- match.arg(method)

    xout <- x

    data.name <- deparse(substitute(x))
    dname <- as.character(substitute(x))
    if(length(dname) > 1) dname <- c(dname[2], dname[length(dname)])
    else dname <- data.name

    u <- threshold

    x <- na.action(x)
    n <- length(x)

    eid <- x > u
    N.u <- sum(eid)

    if(sum(eid)==0) {
	res <- numeric(0)
	attr(res, "theta.msg") <- "extremalindex: No values of x exceed the threshold."
	attr(res, "method") <- "none"
	class(res) <- "extremalindex"
	return(res)
    }

    if(method=="intervals") {

	T.u <- diff((1:n)[eid])

	if(!any(T.u > 2)) {
            hold <- colSums(cbind(T.u, T.u^2))
	    theta.type <- "theta.hat"
	} else {
	    hold <- colSums(cbind(T.u - 1, (T.u - 1) * (T.u - 2)))
	    theta.type <- "theta.tilde"
	}

	theta <- 2 * (hold[1])^2/((N.u - 1) * hold[2])

	if(theta > 1) {
	    theta <- 1
	    theta.msg <- "Extremal index estimate > 1.  Re-set to 1."
	} else theta.msg <- "Valid extremal index estimated"

	K <- ifelse(round(theta * N.u, digits=0) != theta * N.u, ceiling(theta * N.u), theta * N.u)

	o <- order(T.u, na.last=TRUE, decreasing=TRUE)
	oT.u <- T.u[o]
	r <- oT.u[K]

	##
	## Handle ties
	## 
	T.u2 <- c(diff(oT.u), 0)
	ind0 <- T.u2 == 0
	ind1 <- !ind0
	if((K > 1) && (T.u2[K - 1] == 0)) {
	    ind2.0 <- (1:(N.u - 1))[ind0]
	    ind2.1 <- (1:(N.u - 1))[ind1]
	    ind3 <- (ind2.0 < K) & (ind2.0 > max(ind2.1[ind2.1 < K]))
	    K <- min(ind2.0[ind3])
	    r <- oT.u[K]
	}

	res <- c(theta, K, r)
	attr(res, "theta.type") <- theta.type

    } else if(method == "runs") {

	tmp <- decluster(x=x, threshold=threshold, method="runs", r=run.length, ...)
	nc <- length(unique(attributes(tmp)$clusters))
	N.r <- sum(c(tmp) > threshold)
	res <- c(N.r/N.u, nc, run.length)
	attr(res, "cluster") <- attributes(tmp)$cluster

    } else stop("extremalindex: sorry, currently only methods runs and intervals available.")

    names(res) <- c("extremal.index", "number.of.clusters", "run.length")

    attr(res, "data") <- xout
    attr(res, "method") <- method
    attr(res, "data.name") <- dname
    attr(res, "data.call") <- data.name
    attr(res, "call") <- match.call()
    attr(res, "na.action") <- na.action
    attr(res, "threshold") <- u
    class(res) <- "extremalindex"

    return(res)

} # end of 'extremalindex' function.

ci.extremalindex <- function(x, alpha=0.05, R=502, return.samples=FALSE, ...) {

    tmp <- attributes(x)

    r <- x["run.length"]
    nc <- x["number.of.clusters"]

    y <- datagrabber(x)
    n <- length(y)
    u <- tmp$threshold
    eid <- y > u
    Nu <- sum(eid)
    s <- (1:n)[eid]
    Tu <- diff(s)

    if(tmp$method=="runs") cluster <- tmp$cluster
    else if(tmp$method=="intervals") {
	cluster <- rep(1, Nu)
        if(Nu > 1) cluster[2:Nu] <- 1 + cumsum(Tu > r)
    }

    # Set of intercluster times.
    # ict <- c(unlist(aggregate(c(Tu,NA), by=list(cluster), function(x) return(x[1]))$x))
    # ict <- ict[!is.na(ict)]

    # Data by cluster.
    csizes <- tabulate(cluster)
    csizes <- csizes[!is.na(csizes)]
    csizes <- csizes[csizes != 0]
    mcs <- max(csizes)
    Y <- Thold <- numeric(mcs * nc) + NA

    ind <- unlist(c(apply(cbind(csizes, 0:(nc - 1)), 1, function(x, sz) 1:x[1] + x[2] * sz, sz=mcs)))
    Y[ind] <- y[eid]
    Y <- matrix(Y, mcs, nc)
    Thold[ind] <- c(1, Tu)
    Thold <- matrix(Thold, mcs, nc)

    sam <- matrix(NA, R, 3)

    for(i in 1:R) {

	ic <- sample(1:nc, replace=TRUE)
	Tub <- Thold[,ic]
	if(nc > 1) Tub[1,2:nc] <- sample(Thold[1, 2:nc], replace=TRUE)
	# if(nc > 1) Tub[1,2:nc] <- sample(ict, replace=TRUE)

	Tub <- Tub[!is.na(Tub)]
	sb <- cumsum(Tub)

	Yb.star <- Y[,ic]
	Yb <- rep(u, max(sb))
	Yb[sb] <- c(Yb.star[!is.na(Yb.star)])

	sam[i,] <- extremalindex(x=Yb, threshold=u, method=tmp$method, run.length=r)

    } # end of for 'i' loop.

    if(is.null(colnames(sam))) colnames(sam) <- c("extremal.index", "number.of.clusters", "run.length")
    if(tmp$method=="runs") sam <- sam[,-3]

    if(return.samples) return(sam)

    bds <- t(apply(sam, 2, quantile, probs=c(alpha/2, 1 - alpha/2)))

    conf.level <- paste((1 - alpha) * 100, "%", sep="")
    if(tmp$method=="runs") out <- cbind(bds[,1], x[1:2], bds[,2])
    else out <- cbind(bds[,1], x, bds[,2])
    colnames(out) <- c(paste(conf.level, " lower CI", sep=""), "Estimate", paste(conf.level, " upper CI", sep=""))

    attr(out, "data.name") <- tmp$call
    attr(out, "method") <- "Ferro-Segers Bootstrap"
    attr(out, "ei.method") <- tmp$method
    if(tmp$method=="runs") {
	names(r) <- "run.length"
	attr(out, "parameter") <- r
    }
    attr(out, "conf.level") <- conf.level
    class(out) <- "ci"
    return(out)
} # end of 'ci.extremalindex' function.

print.extremalindex <- function(x, ...) {
    tmp <- attributes(x)
    if(tmp$method=="intervals") {
	cat("\n", "Intervals Method Estimator for the Extremal Index\n")
	print(tmp$theta.msg)
	if(tmp$theta.type=="theta.hat") cat("\n", "theta.hat used because there are no inter-exceedance times > 2.\n") 
	else if(tmp$theta.type=="theta.tilde") cat("\n", "theta.tilde used because there exist inter-exceedance times > 2.\n")
    } else if(tmp$method=="runs") {
	cat("\n", "Runs Estimator for the Extremal Index\n")
    } else if(tmp$method=="none") {
	print(tmp$theta.msg)
    } else stop("print: sorry, currently no print method for anything but intervals or runs estimation methods.")
    y <- x
    attributes(y) <- NULL
    names(y) <- tmp$names
    print(y)
    invisible()
} # end of 'print.extremalindex' function.

decluster <- function(x, threshold, ...) {
    UseMethod("decluster", x)
} # end of 'decluster' function.

# decluster.formula <- function(x, threshold, ..., method=c("runs","intervals"), clusterfun="max") {
#    method <- tolower(method)
#    method <- match.arg(method)
#    class(x) <- c(method, class(x))
#    UseMethod("decluster", x)
# }

decluster.default <- function(x, threshold, ..., method=c("runs","intervals"), clusterfun="max") {

    method <- tolower(method)
    method <- match.arg(method)
    class(x) <- c(method, class(x))
    UseMethod("decluster", x)

} # end of 'decluster.default' function.

decluster.data.frame <- function(x, threshold, ..., which.cols, method=c("runs","intervals"), clusterfun="max") {

    if(missing(which.cols)) stop("decluster: Must provide a which.cols argument for data frames.")
    if(!is.element(length(which.cols), 1:2)) stop("decluster: which.cols must have length one or two.")

    y <- x[,which.cols[1]]
    if(length(which.cols) == 2) groups <- x[,which.cols[2]]
    else groups <- NULL

    out <- decluster(x=y, threshold=threshold, method=method, clusterfun=clusterfun, groups=groups, ...)
    attr(out, "data.name") <- c(deparse(substitute(x)), which.cols)
    if(is.null(groups)) attr(out, "groups") <- "NULL"

    return(out)

} # end of 'decluster.data.frame' function.

decluster.runs <- function(x, threshold, ..., data, r=1, clusterfun="max", groups=NULL, replace.with, na.action=na.fail) {

    clusterfun <- match.fun(clusterfun)
    if(missing(replace.with)) replace.with <- threshold

    xout <- x

    dname <- deparse(substitute(x))
    thname <- deparse(substitute(threshold))

#     if(!missing(data)) {
# 	# if(missing(data)) stop("decluster: must supply data argument if x is a formula.")
# 	datc <- colnames(data)
# 	if(is.element(dname, datc)) {
# 	    id <- dname == datc
# 	    x.fun <- ifelse(dname == "substitute(x)", deparse(x), dname)
# 	    x <- c(data[,id])
#             x.fun <- formula(paste(x.fun, "~ 1"))
# 	} else if(is.formula(x)) {
# 	    x.fun <- x
# 	    x.str <- unlist(strsplit(as.character(x.fun), split="~"))
# 	    x.str <- x.str[x.str != ""]
# 	    ns <- length(x.str)
# 	    if(x.str[ns] != 1) stop("decluster: invalid x argument.")
# 	    x <- eval(parse(text=x.str[1]), envir=data)
#         } else x.fun <- NULL
# 
# 	if(is.element(thname, datc)) {
# 	    id <- thname == datc
# 	    th.fun <- ifelse(thname == "substitute(theshold)", deparse(threshold), thname)
# 	    threshold <- c(data[,id])
# 	    th.fun <- formula(paste(th.fun, "~ 1"))
# 	} else if(is.formula(threshold)) {
# 	    th.fun <- threshold
# 	    th.str <- unlist(strsplit(as.character(th.fun), split="~"))
# 	    th.str <- th.str[th.str != ""]
# 	    ns <- length(th.str)
# 	    if(th.str[ns] != 1) stop("decluster: invalid threshold argument.")
# 	    threshold <- eval(parse(text=th.str[1]), envir=data)
#         }
#     } else {
# 	if(is.formula(x) || is.formula(threshold)) stop("decluster: data argument missing.")
# 	x.fun <- th.fun <- NULL
#     }

    if(length(threshold==1)) x <- na.action(x)
    else {
	look <- na.action(cbind(x, threshold))
	x <- look[,1]
	threshold <- look[,2]
    }
    n <- length(x)

    eid <- x > threshold
    xu <- x[eid]

    Nu <- sum(eid)
    s <- (1:n)[eid]
    cluster <- rep(1, Nu)

    if(is.null(groups)) {

	res <- numeric(n) + replace.with
	res[!eid] <- x[!eid]

	T.u <- diff(s)
	if(Nu > 1) {

	    cluster[2:Nu] <- 1 + cumsum(T.u > r)

	    # For some unknown reason, the next command didn't work, but the subsequent one does.
	    # cres <- c(aggregate(xu, by=list(cluster), clusterfun, ...)$x)[cluster]
	    cres <- c(aggregate(xu, by=list(cluster), clusterfun)$x)[cluster]
	    cf.id <- xu == cres

	    posfun <- function(ind) {
	        m <- length(ind)
	        res <- logical(m)
	        if(any(ind)) {
	        	ind2 <- min((1:m)[ind])
	        	res[ind2] <- TRUE
	        } else res[1] <- TRUE
	        return(res)
	    } # end of internal 'posfun' function.
	    pos <- ((1:n)[eid])[unlist(c(t(aggregate(cf.id, by=list(cluster), posfun)$x)))]
	    res[pos] <- x[pos]
	} else if(any(eid)) res[eid] <- x[eid]
	
    } else {

	b <- groups
	if(length(threshold) == 1) {

	    res <- aggregate(x, by=list(b), FUN=decluster, threshold=threshold, method="runs", r=r, clusterfun=clusterfun, groups=NULL)
	    res <- unlist(c(t(res$x)))

	} else {

	    dcfun <- function(dat, r, clusterfun) {
		xdat <- dat[,1]
		th <- dat[,2]
		return(decluster(x=xdat, threshold=th, method="runs", r=r, clusterfun=clusterfun))
	    } # end of internal 'dcfun' function.

	    hold <- list()
	    ub <- unique(b)
	    nb <- length(ub)
	    for(i in 1:nb) hold[[i]] <- cbind(x[b==ub[i]], threshold[b==ub[i]])
	    res <- lapply(hold, FUN=dcfun, r=r, clusterfun=cluster)
	    res <- c(unlist(res))
	} # end of if else 'threshold' is constant or not stmts.

	T.u <- diff(s)
	k <- c(aggregate(groups, by=list(groups), length)$x)
	bl <- rep(1:length(unique(groups)), k)[eid] - 1
        if(Nu > 1) cluster[2:Nu] <- 1 + cumsum(T.u > r) + bl[2:Nu]

    } # end of if else 'groups' are NULL stmts.

    attr(res, "call") <- match.call()
    # attr(res, "x.fun") <- x.fun
    # attr(res, "th.fun") <- th.fun
    attr(res, "data") <- xout
    attr(res, "data.name") <- dname
    attr(res, "decluster.function") <- clusterfun
    attr(res, "method") <- "runs"
    attr(res, "threshold") <- threshold
    attr(res, "groups.name") <- deparse(substitute(groups))
    attr(res, "groups") <- groups
    attr(res, "run.length") <- r
    attr(res, "na.action") <- na.action
    attr(res, "clusters") <- cluster
    class(res) <- "declustered" 

    return(res)

} # end of 'decluster.runs' function.

decluster.intervals <- function(x, threshold, ..., clusterfun="max", groups=NULL, replace.with, na.action=na.fail) {

    xout <- x

#     if(is.formula(x)) {
#         if(missing(data)) stop("decluster: must supply data argument if x is a formula.")
#         x.fun <- x
# 	x.fun <- ifelse(deparse(substitute(x)) == "substitute(x)", deparse(x), deparse(substitute(x)))
#         # x.fun <- formula(paste(x.fun, "~ 1"))
#         x <- model.response(model.matrix(x.fun, data=data))
#     } else x.fun <- NULL
# 
#     if(is.formula(threshold)) {
#         if(missing(data)) stop("decluster: must supply data argument if threshold is a formula.")
#         th.fun <- threshold
# 	th.fun <- ifelse(deparse(substitute(threshold)) == "substitute(threshold)", deparse(threshold), deparse(substitute(threshold)))
#         # th.fun <- formula(paste(th.fun, "~ 1"))
#         threshold <- model.response(model.matrix(th.fun, data=data))
#     } else th.fun <- NULL

    th <- extremalindex(x=x, threshold=threshold, na.action=na.action)
    if(th[1] >= 1) r <- 0
    else r <- th[3]

    if(missing(replace.with)) tmp <- decluster(x=x, threshold=threshold, clusterfun=clusterfun, groups=groups, na.action=na.action)
    else tmp <- decluster(x=x, threshold=threshold, clusterfun=clusterfun, groups=groups, replace.with=replace.with, na.action=na.action)

    hold <- attributes(tmp)
    res <- c(tmp)

    attr(res, "call") <- match.call()
    # attr(res, "x.fun") <- x.fun
    # attr(res, "th.fun") <- th.fun
    attr(res, "data") <- xout
    attr(res, "data.name") <- deparse(substitute(x))
    attr(res, "decluster.function") <- clusterfun
    attr(res, "method") <- "intervals"
    attr(res, "threshold") <- threshold
    attr(res, "groups.name") <- deparse(substitute(groups))
    attr(res, "groups") <- groups
    attr(res, "run.length") <- r
    attr(res, "na.action") <- na.action
    attr(res, "clusters") <- hold$clusters
    class(res) <- "declustered"

    return(res)

} # end of 'decluster.intervals' function.


plot.declustered <- function(x, which.plot=c("scatter", "atdf"), qu=0.85, xlab=NULL, ylab=NULL, main=NULL, col="gray", ...) {

    which.plot <- tolower(which.plot)
    which.plot <- match.arg(which.plot)

    tmp <- attributes(x)
     if(length(tmp$data.name) == 1) {
 	# x0 <- try(get(tmp$data.name), silent=TRUE)
 	# if(class(x0) == "try-error") x0 <- NULL
	x0 <- tmp$data
     } else x0 <- datagrabber(x)

    if(!is.null(x0)) {
        if(is.data.frame(x0) || is.matrix(x0)) {
	    bl <- x0[,2]
	    x0 <- x0[,1]
        } else bl <- NULL
    } else bl <- NULL

    x <- c(x)
    # if(tmp$data.name == "x") warning("plot.declustered: attempt to plot original data failed because original data has same name (x) as declustered argument.")
    # x0 <- tmp$na.action(x0)
    if(!is.null(x0)) n <- length(x0)
    else n <- length(x)
    m <- length(x)
 
    if(is.null(xlab)) xlab <- "index"
    if(is.null(ylab)) {
	ylab <- tmp$data.name
	if(length(ylab) == 3) ylab <- paste(ylab[1], ": ", ylab[2], " grouped by ", ylab[3], sep="")
	else if(length(ylab) == 2) ylab <- paste(ylab[1], ": ", ylab[2], sep="")
    }
    if(is.null(main) && which.plot == "scatter") main <- deparse(tmp$call)
    else main <- ""

    # if(which.plot == "primary") par(mfrow=c(5,1), oma=c(0,0,2,0), mar=c(5.1, 4, 0.1, 2.1))
    # else
    if(which.plot == "atdf") {
	if(!is.null(x0)) par(mfrow=c(2,2), oma=c(0,0,2,0), mar=c(4.1, 4.1, 2.1, 2.1))
	else par(mfrow=c(2,1), oma=c(0,0,2,0), mar=c(4.1, 4.1, 2.1, 2.1))
    }

    if(is.element(which.plot, c("primary", "scatter"))) {
        if(!is.null(x0)) plot(x0, xlab=xlab, ylab=ylab, main=main, col=col, ...)
	else plot(x, xlab=xlab, ylab=ylab, main=main, col=col, ...)
        u <- tmp$threshold
        if(length(u) == 1 || all(u == u[1])) {
	    u <- u[1]
	    abline(h=u, lty=2)
        } else if(length(u) == n) lines(1:n, u, lty=2)
        else lines(1:length(u), u, lty=2)
        if(tmp$groups.name != "NULL") {
	    if(is.null(bl)) bl <- tmp$groups
	    if(class(bl) != "try-error") # warning("plot.declustered: unable to get groups for clusters.")
	    {
		m2 <- length(unique(bl))
	        m3 <- length(bl)
	        k <- c(aggregate(bl, by=list(bl), length)$x)
	        bl <- rep(1:m2, k)
	        bl2 <- (1:m3)[diff(bl) == 1]
	        abline(v=bl2, col="blue")
	    }
        }
        points(1:m, x)
    }

    if(is.element(which.plot, c("primary", "atdf"))) {
	if(!is.null(x0)) hold0 <- atdf(x0, u=qu, plot=FALSE, na.action=tmp$na.action)
	hold1 <- atdf(x, u=qu, plot=FALSE, na.action=tmp$na.action)
	if(which.plot=="atdf") {
	    if(!is.null(x0)) plot(hold0, type="rho", ylab="rho (original data)", xlab="", main=paste(ylab[1], " (original)\n", "auto-tail dependence function", sep=""))
	    plot(hold1, type="rho", ylab="rho (declustered data)", xlab="", main=paste(ylab[1], " (declustered)\n", "auto-tail dependence function", sep=""))
 	} else {
	    if(!is.null(x0)) plot(hold0, type="rho", ylab="rho (original data)", xlab="", main="")
            plot(hold1, type="rho", ylab="rho (declustered data)", xlab="", main="")
	}
	if(which.plot == "primary" && !is.null(x0)) plot(hold0, type="rhobar", ylab="rhobar (original data)", xlab="", main="")
	else if(!is.null(x0)) plot(hold0, type="rhobar", ylab="rhobar (original data)", main="")
	plot(hold1, type="rhobar", ylab="rhobar (declustered data)", main="")
    }

    if(is.element(which.plot, c("primary","atdf"))) mtext(tmp$call, line=0.5, outer=TRUE)
    invisible()
} # end of 'plot.declustered' function.

print.declustered <- function(x, ...) {
    tmp <- attributes(x)
    cat("\n", tmp$data.name, " declustered via", tmp$method, " declustering.\n")
    hold <- extremalindex(x, threshold=tmp$threshold)
    cat("\n", "Estimated extremal index (intervals estimate) = ", hold[1], "\n")
    cat("\n", "Number of clusters = ", length(unique(tmp$clusters)), "\n")
    cat("\n", "Run length = ", tmp$run.length, "\n")
    invisible()
} # end of 'print.declustered' function.
