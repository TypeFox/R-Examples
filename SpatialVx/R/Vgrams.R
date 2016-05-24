griddedVgram <- function(object, zero.in=TRUE, zero.out=TRUE, time.point=1, model=1, ...) {

    a <- attributes(object)
    out <- list()
    if(!is.null(a$names)) a$names <- NULL
    attributes(out) <- a

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(object, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(object, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(object, model=model)
    else dat <- datagrabber(object)

    X <- dat$X
    Y <- dat$Xhat
    ## End: Get the data sets

    out$zero.in <- zero.in
    out$zero.out <- zero.out

    out$Vx.vgram.matrix <- list()
    out$Fcst.vgram.matrix <- list()

    if(zero.in) {
	out$Vx.vgram.matrix[[1]] <- vgram.matrix(dat=X, ...)
	out$Fcst.vgram.matrix[[1]] <- vgram.matrix(dat=Y, ...)
    } else {
	out$Vx.vgram.matrix[[1]] <- NULL
	out$Fcst.vgram.matrix[[1]] <- NULL
    }

    if(zero.out) {
	out$Vx.vgram.matrix[[2]] <- variogram.matrix(dat=X, zero.out=TRUE, ...)
	out$Fcst.vgram.matrix[[2]] <- variogram.matrix(dat=Y, zero.out=TRUE, ...)
    } else {
	out$Vx.vgram.matrix[[2]] <- NULL
	out$Fcst.vgram.matrix[[2]] <- NULL
    }

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[1:2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    attr(out, "data.name") <- c(vxname, dn[model.num])

    class(out) <- "griddedVgram"
    return(out)
} # end of 'griddedVgram' function.

plot.griddedVgram <- function(x, ..., set.pw=FALSE) {

    a <- attributes(x)

    if(length(a$data.name) == 3) {
	mainX <- a$data.name[2]
        mainY <- a$data.name[3]
    } else {
        mainX <- a$data.name[1]
        mainY <- a$data.name[2]
    }

   if(x$zero.in) {
      vgX <- x$Vx.vgram.matrix[[1]]
      vgY <- x$Fcst.vgram.matrix[[1]]
   }

   if(x$zero.out) {
      vgX.zero <- x$Vx.vgram.matrix[[2]]
      vgY.zero <- x$Fcst.vgram.matrix[[2]]
   }

    if(!is.logical(set.pw) && !is.numeric(set.pw)) stop("plot.griddedVgram: invalid set.pw argument.")
    else if(!is.logical(set.pw) && is.numeric(set.pw)) {

	if(length(set.pw) != 2) stop("plot.griddedVgram: invalid set.pw argument.")
	par(mfrow=set.pw, oma=c(0,0,2,0))

    } else if(set.pw) {

	if(((x$zero.in) && !(x$zero.out)) || (!(x$zero.in) & (x$zero.out))) {
	    par(mfrow=c(2,2), mar=rep(4.1,4), oma=c(0,0,2,0))
	} else if(x$zero.in && x$zero.out) par(mfrow=c(4,2), mar=rep(4.1,4), oma=c(0,0,2,0))

    } else par(oma=c(0,0,2,0))

   if(x$zero.in) {

      plot(vgX$d, vgX$vgram, xlab="separation distance", ylab="variogram", main=mainX, col="darkblue")
      points(vgX$d.full, vgX$vgram.full, pch=".", cex=1.25, col="darkblue")

      plot(vgY$d, vgY$vgram, xlab="separation distance", ylab="variogram", main=mainY, col="darkblue") 
      points(vgY$d.full, vgY$vgram.full, pch=".", cex=1.25, col="darkblue")

      class(vgX) <- "vgram.matrix"
      class(vgY) <- "vgram.matrix"

      zl <- range(c(c(vgX$vgram.full), c(vgY$vgram.full)), finite=TRUE)
      plot(vgX, xlab="x separations", ylab="y separations", zlim=zl, ...)
      plot(vgY, xlab="x separations", ylab="y separations", zlim=zl, ...)

   }
   if(x$zero.out) {
      plot(vgX.zero$d, vgX.zero$vgram, xlab="separation distance", ylab="variogram", main=paste(mainX, " (non-zeros only)", sep=""), col="darkblue")
      points(vgX.zero$d.full, vgX.zero$vgram.full, pch=".", cex=1.25, col="darkblue")
      plot(vgY.zero$d, vgY.zero$vgram, xlab="separation distance", ylab="variogram", main=paste(mainY, " (non-zeros only)", sep=""), col="darkblue")
      points(vgY.zero$d.full, vgY.zero$vgram.full, pch=".", cex=1.25, col="darkblue")
      class(vgX.zero) <- "vgram.matrix"
      class(vgY.zero) <- "vgram.matrix"
      zl <- range(c(c(vgX.zero$vgram.full), c(vgY.zero$vgram.full)), finite=TRUE)
      plot(vgX.zero, xlab="x separations", ylab="y separations", zlim=zl, ...)
      plot(vgY.zero, xlab="x separations", ylab="y separations", zlim=zl, ...)
   }

    if(!is.null(a$msg)) {
	title("")
	mtext(a$msg, line=0.05, outer=TRUE)
    }

   invisible()
} # end of 'plot.griddedVgram' function.

corrskill <- function(x,y,...) {
   good <- !is.na(x) & !is.na(y)
   n <- sum(good, na.rm=TRUE)
   s1 <- sd(c(x), ...)
   s2 <- sd(c(y), ...)
   m1 <- sum(colSums(x,na.rm=TRUE),na.rm=TRUE)/n
   m2 <- sum(colSums(y,na.rm=TRUE),na.rm=TRUE)/n
   out <- ((n/(n-1))/(s1*s2))*(x-m1)*(y-m2)
   return(out)
} # end of 'corrskill' function.

sqerrloss <- function(x,y,...) return((x-y)^2)
abserrloss <- function(x,y,...) return(abs(x-y))
distmaploss <- function(x,y, threshold=0, const=Inf, ...) {
   if(length(threshold)==1) threshold <- rep(threshold,2)
   xdim <- dim(x)
   if(is.finite(const)) {
	x <- cbind(matrix(0, xdim[1], const), x, matrix(0, xdim[1], const))
	x <- rbind(matrix(0, const, xdim[2]+2*const), x, matrix(0, const, xdim[2] + 2*const))
	y <- cbind(matrix(0, xdim[1], const), y, matrix(0, xdim[1], const))
        y <- rbind(matrix(0, const, xdim[2]+2*const), y, matrix(0, const, xdim[2] + 2*const))
   }
   x <- im(x)
   y <- im(y)
   x <- solutionset(x >= threshold[2])
   y <- solutionset(y >= threshold[1])
   bb <- boundingbox(as.rectangle(x), as.rectangle(y))
   x <- rebound(x, bb)
   y <- rebound(y, bb)
   dy <- distmap(y, ...)$v
   dx <- distmap(x, ...)$v
   if(is.finite(const)) {
	dy[dy>const] <- const
	dx[dx>const] <- const
	dy <- dy[(const+1):(xdim[1]+const),(const+1):(xdim[2]+const)]
	dx <- dx[(const+1):(xdim[1]+const),(const+1):(xdim[2]+const)]
   }
   return(abs(dx - dy))
} # end of 'distmaploss' function.

lossdiff <- function(x, ...) {

    UseMethod("lossdiff", x)

} # end of 'lossdiff' function.

lossdiff.SpatialVx <- function(x, ..., time.point = 1, model = c(1,2), threshold = NULL,
    lossfun = "corrskill", zero.out = FALSE) {

    theCall <- match.call()

    a <- attributes(x)

    # if( missing( maxrad ) ) maxrad <- floor( max( a$xdim ) / 2 )

    # if( is.null( maxrad ) || is.na( maxrad ) || !is.finite( maxrad ) || !is.numeric( maxrad ) ) maxrad <- 20

    ## Begin: Get the data sets
    if(!missing(time.point)) dat <- datagrabber(x, time.point=time.point, model=model[1])
    else dat <- datagrabber(x, model=model[1])
   
    X <- dat$X
    Xhat <- dat$Xhat

    if(!missing(time.point)) dat <- datagrabber(x, time.point=time.point, model=model[2])
    else dat <- datagrabber(x, model=model[2])

    Xhat2 <- dat$Xhat
    ## End: Get the data sets

    out <- lossdiff.default(x = X, ..., xhat1 = Xhat, xhat2 = Xhat2, threshold = threshold, lossfun = lossfun, 
			    loc = a$loc, zero.out = zero.out)

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    out$data.name <- c(vxname, dn[model.num])
    out$call <- theCall

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model

    attr(out, "msg") <- a$msg
    attr(out, "projection") <- a$projection
    attr(out, "map") <- a$map
    attr(out, "loc") <- a$loc
    attr(out, "loc.byrow") <- a$loc.byrow
    attr(out, "xdim") <- a$xdim

    return(out)

} # end of 'lossdiff.SpatialVx' function.

lossdiff.default <- function(x, ..., xhat1, xhat2, threshold = NULL, lossfun = "corrskill", loc = NULL, zero.out = FALSE) {

    theCall <- match.call()

    out <- list()
    attributes(out) <- atmp <- list(...)
    if(is.null(atmp$loc.byrow)) attr(out, "loc.byrow") <- FALSE

    data.name <- c(as.character(substitute(x)),as.character(substitute(xhat1)),as.character(substitute(xhat2)))
    names(data.name) <- c("verification","model1", "model2")
    out$data.name <- data.name
    xdim <- dim(x)

    out$maxdim <- floor( sqrt(xdim[ 1 ]^2 + xdim[2]^2 ) )
    out$lossfun <- lossfun
    out$lossfun.args <- list(...)

    if( !is.null( threshold ) ) {

	numthresh <- length( threshold )
	if( numthresh == 1 ) threshold <- rep( threshold, 3 )
	else if( numthresh > 3 ) stop("lossdiff: invalid threshold argument.  Must be of length one, two or three.")
	else if( numthresh == 2 ) threshold <- c( threshold[ 1 ], rep( threshold[ 2], 2 ) )

	x[ x < threshold[ 1 ] ] <- 0
	xhat1[ xhat1 < threshold[ 2 ] ] <- 0
	xhat2[ xhat2 < threshold[ 3 ] ] <- 0

	out$threshold <- threshold

    } # end of if 'threshold' stmt.
 
    g1 <- do.call(lossfun, args=c(list(x=x,y=xhat1),list(...)))
    g2 <- do.call(lossfun, args=c(list(x=x,y=xhat2),list(...)))
    d <- matrix(g1-g2, xdim[1], xdim[2])

    if(zero.out) {

	# zeros <- d == 0
	zeros <- (x==0) & (xhat1==0) & (xhat2==0)
	beta <- mean(!zeros, na.rm=TRUE)
	out$zeros <- zeros
	out$beta <- beta

    }

#     if(trend=="ols") {
 
 	if(is.null(loc)) loc <- cbind(rep(1:xdim[1],xdim[2]), rep(1:xdim[2],each=xdim[1]))
 	# else out$loc <- character(substitute(loc))
 	if( zero.out ) dat <- data.frame(y = c(d)[ !zeros ], x1 = loc[ !zeros, 1], x2 = loc[ !zeros, 2])
	else dat <- data.frame(y = c(d), x1 = loc[, 1], x2 = loc[, 2] )
 	fit <- try( lm(y~x1+x2, data=dat), silent = TRUE )
 	# tr <- matrix(predict(fit),xdim[1],xdim[2])
 	# if(zero.out) tr[zeros] <- 0
 	# d <- d - tr
 	out$trend.fit <- fit
# 
#     } else if(is.numeric(trend)) {
#
#	if(zeros) warning("lossdiff: zero.out is TRUE, but trend provided.  Not setting original zeros back after removing trend.")
# 	d <- d - trend
# 	if( zero.out ) d[ zero.out ] <- 0
# 	out$trend <- trend
#
#   }

    out$call <- theCall
    out$d <- d
    class(out) <- "lossdiff"

    return(out)

} # end of 'lossdiff.default' function.

empiricalVG.lossdiff <- function( x, trend = 0, maxrad, dx = 1, dy = 1 ) {

    theCall <- match.call()

    a <- attributes( x )

    out <- x

    if( missing( maxrad ) ) maxrad <- floor( max( a$xdim ) / 2 )

    if( is.null( maxrad ) || is.na( maxrad ) || !is.finite( maxrad ) || !is.numeric( maxrad ) ) maxrad <- 20

    out$trend <- trend
    out$vgram.args <- list( maxrad = maxrad, dx = dx, dy = dy )

    if( !is.null( x$zeros ) ) zero.out <- TRUE
    else zero.out <- FALSE

    d <- x$d - trend

    if( zero.out ) d[ x$zeros ] <- 0

    out$loss.differential.detrended <- d

    if( !zero.out ) vg <- vgram.matrix(dat = d, R = maxrad, dx = dx, dy = dy)
    else vg <- variogram.matrix(dat = d, R = maxrad, dx = dx, dy = dy, zero.out = zero.out)

    out$lossdiff.vgram <- vg

    out$call <- list( out$call, theCall )

    return( out )

} # end of 'empiricalVG.lossdiff' function.

flossdiff <- function(object, vgmodel = "expvg", ... ) {

    out <- object
    a <- attributes( object )
    vg <- object$lossdiff.vgram

    p <- c( sqrt( vg$vgram[ 1 ] ), ifelse( vg$vgram[ 2 ] > 0, -(vg$vgram[ 2 ] - vg$vgram[ 1 ]), -(0 - vg$vgram[ 1 ]) ) )

    fit <- try( nlminb( p, ORSS, vg = vg, model = vgmodel, ..., lower = c(0, 0), upper = c(Inf, Inf) ) )

    hold <- list( model = vgmodel, p = fit$par, objective = fit$objective, maxdim = object$maxdim,
		convergence = fit$convergence, message = fit$message, iterations = fit$iterations,
		evaluations = fit$evaluations )

    class( hold ) <- paste("flossdiff", vgmodel, sep = ".")

    out$vgmodel <- hold

    # class( out ) <- c("lossdiff", "flossdiff")

    return(out)

} # end of 'flossdiff' function.

expvg <- function( p, vg, ... ) {

    return( p[ 1 ] * ( 1 - exp( - vg$d * p[ 2 ] ) ) )

} # end of 'expvg' function.

predict.flossdiff.expvg <- function( object, newdata, ... ) {

    if( missing( newdata ) ) {

	h <- seq(0, object$maxdim - 1, by = 1 )

    } else h <- newdata

    p <- object$p

    return( p[ 1 ] * ( 1 - exp( - h * p[ 2 ] ) ) )

} # end of 'predict.expvg' function.

print.flossdiff.expvg <- function( x, ... ) {

    cat("\n\nExponential variogram fit\n")

    look <- x$p
    names( look ) <- c("scale", "range")

    print( look )

    invisible()

} # end of 'print.flossdiff.expvg' function.

ORSS <- function( p, vg, model = "expvg", ... ) {

    fit <- do.call( model, c( list( p = p, vg = vg ), list( ... ) ) )

    out <- sum( ( fit - vg$vgram )^2, na.rm = TRUE )

    return( out )

} # end of 'ORSS' function.

# TO DO: make a WRSS function too for weighted RSS (maybe not necessary).

plot.lossdiff <- function(x, ..., icol=c("gray", tim.colors(64))) {

    tmp <- attributes(x)
    loc.byrow <- tmp$loc.byrow

    if(is.null(tmp$msg)) msg <- paste("\n", x$lossfun, ": ", x$data.name[2], " vs ", x$data.name[3], " (", x$data.name[1], ")", sep="")
    else msg <- paste(tmp$msg, "\n", x$lossfun, ": ", x$data.name[2], " vs ", x$data.name[3], " (", x$data.name[1], ")", sep="")
    
    par(mfrow=c(2,2), mar=rep(4.1,4), oma=c(0,0,2,0))

    if(is.null(tmp$projection)) proj <- FALSE
    else proj <- tmp$projection

    if(is.null(tmp$map)) domap <- FALSE
    else domap <- tmp$map

    if(is.null(tmp$xdim)) proj <- domap <- FALSE
    else xd <- tmp$xdim

    if(proj) {
	loc <- list(x=matrix(tmp$loc[,1], xd[1], xd[2], byrow=loc.byrow),
		    y=matrix(tmp$loc[,2], xd[1], xd[2], byrow=loc.byrow))
    }

    if(domap) {
	locr <- apply(tmp$loc, 2, range, finite=TRUE)
	ax <- list(x=pretty(round(tmp$loc[,1], digits=2)), y=pretty(round(tmp$loc[,2], digits=2)))
    }

    if(is.null(x$zeros)) Im <- x$d
    else {
        Im <- x$d
        Im[x$zeros] <- NA
    }

    if(domap) {
	map(xlim=locr[,1], ylim=locr[,2], type="n")
	axis(1, at=ax$x, labels=ax$x)
	axis(2, at=ax$y, labels=ax$y)
        if(proj) image.plot(loc$x, loc$y, Im, col=icol, add=TRUE, ...)
        else image.plot(Im, col=icol, add=TRUE, ...)
        map(add=TRUE, lwd=1.5)
        map(add=TRUE, database="state")
    } else {
	if(proj) image.plot(loc$x, loc$y, Im, col=icol, ...)
        else image.plot(Im, col=icol, ...)
    }
    title("")
    title("Loss Differential Field")

    hist(x$d, breaks="FD", xlab="Mean Loss Differential", col="darkblue", freq=FALSE,
	main="Histogram of\nLoss Differential Field")
    a <- x$lossdiff.vgram
    plot(a$d, a$vgram, col="darkblue", xlab="separation distance", ylab="variogram")
    if(!is.null(x$vgmodel)) {
	lines(a$d, predict(x$vgmodel, newdata = a$d), col="darkorange", lwd=1.5)
   	legend("bottomright", legend=c("Empirical", "Model"), pch=c("o", ""), col=c("darkblue","darkorange"), lty=c(0,1), lwd=1.5, bty="n")
    }
    plot.vgram.matrix(a, main="variogram by direction")

    title("")
    mtext(msg, line=0.05, outer=TRUE)

    invisible()

} # end of 'plot.lossdiff' function.

print.lossdiff <- function( x, ... ) {

    print( x$call )

    print( x$data.name )

    if( !is.null( x$trend ) ) print( summary( x$trend ) )

    print( x$lossfun )

    if( !is.null( x$vgmodel ) ) {

	cat("\n", "Fitted variogram model:\n")

	print( x$vgmodel )

    } # end of if model has been fit function.

    invisible()

} # end of 'print.lossdiff' function.

summary.lossdiff <- function(object, ...) {

   out <- object
   out$summary.call <- match.call()
   msg <- paste(object$lossfun, ": ", object$data.name[2], " vs ", object$data.name[3], " (against verification: ", object$data.name[1], ")", sep="")
   print(msg)
   cat("Fitted trend information (note: not used by subsequent functions, information only):\n\n")
   print( summary( object$trend.fit ) )
   d <- object$d
   if(zero.out <- !is.null(object$beta)) {
	cat("\n", "Estimate of beta present, so calculating Dbar and test statistic over non-zero entries only.\n")
	cat("\n", "Frequency of non-zero loss differential (beta) is: ", object$beta, "\n")
	d[d==0] <- NA
   }
   good <- !is.na(d)
   n <- length(d[good])
   cat("\n", "number of non-zero loss differential points is: ", n, "\n")
   cat("\n", "Mean Loss Differential: \n")
   Dbar <- sum(colSums(d,na.rm=TRUE),na.rm=TRUE)/n
   out$Dbar <- Dbar
   cat("\n", Dbar, "\n", "\n")
   a <- object$lossdiff.vgram
   cat("\n", "Summary of empirical variogram values:\n")
   print(stats(a))
   if(!is.null(object$vgmodel)) {
      # co <- coef(object$vgmodel)
      # sig2 <- co[1]^2
      # r <- co[2]
      # b <- sig2*exp(-a$d.full/r)
      # if(zero.out) denom <- sqrt(mean(object$beta^2*b,na.rm=TRUE))
      # else 
      b <- predict(object$vgmodel)
      denom <- sqrt(mean(b,na.rm=TRUE))
      SV <- Dbar/denom
      out$test.statistic <- SV
      cat("Test Statistic for null hypothesis of equal predictive ability on average\n")
      print(SV)
      pval <- c( 2 * pnorm( -abs(SV) ), pnorm( SV ), pnorm( SV, lower.tail=FALSE ) )
      names(pval) <- c("two.sided", "less", "greater")
      cat("p-value for two-sided alternative hypothesis is: ", pval[1], "\n")
      cat("p-value for (one-sided) alternative hypothesis that mu(D) < 0 is: ", pval[2], "\n")
      cat("p-value for (one-sided) alternative hypothesis that mu(D) > 0 is: ", pval[3], "\n")
      out$p.value <- pval
   }

   invisible(out)

} # end of 'summary.lossdiff' function.

variogram.matrix <- function (dat, R = 5, dx = 1, dy = 1, zero.out = FALSE) 
{
    SI <- function(ntemp, delta) {
        n1 <- 1:ntemp
        n2 <- n1 + delta
        good <- (n2 >= 1) & (n2 <= ntemp)
        cbind(n1[good], n2[good])
    }
    if(zero.out) dat[dat==0] <- NA
    N <- ncol(dat)
    M <- nrow(dat)
    m <- min(c(round(R/dx), M))
    n <- min(c(round(R/dy), N))
    ind <- rbind(as.matrix(expand.grid(0, 1:n)), as.matrix(expand.grid(1:m, 
        0)), as.matrix(expand.grid(c(-(m:1), 1:m), 1:n)))
    d <- sqrt((dx * ind[, 1])^2 + (dy * ind[, 2])^2)
    good <- (d > 0) & (d <= R)
    ind <- ind[good, ]
    d <- d[good]
    ind <- ind[order(d), ]
    d <- sort(d)
    nbin <- nrow(ind)
    holdVG <- rep(NA, nbin)
    holdN <- rep(NA, nbin)
    for (k in 1:nbin) {
        MM <- SI(M, ind[k, 1])
        NN <- SI(N, ind[k, 2])
	numNA <- sum(is.na(dat[MM[,1],NN[,1]]) | is.na(dat[MM[,2],NN[,2]]),na.rm=TRUE)
        holdN[k] <- length(MM) * length(NN) - numNA
        BigDiff <- (dat[MM[, 1], NN[, 1]] - dat[MM[, 2], NN[,2]])
        holdVG[k] <- mean(0.5 * (BigDiff)^2, na.rm=TRUE)
    }
    top <- tapply(holdVG * holdN, d, FUN = "sum")
    bottom <- tapply(holdN, d, FUN = "sum")
    dcollapsed <- as.numeric(names(bottom))
    vgram <- top/bottom
    dimnames(vgram) <- NULL
    out <- list(vgram = vgram, d = dcollapsed, ind = ind, d.full = d, 
        vgram.full = holdVG, N = holdN, dx = dx, dy = dy)
    class(out) <- "vgram.matrix"
    return(out)
} # end of 'variogram.matrix' function.

structurogram.matrix <- function(dat, q=2, R=5, dx=1, dy=1, zero.out=FALSE) {
   SI <- function(ntemp, delta) {
        n1 <- 1:ntemp
        n2 <- n1 + delta
        good <- (n2 >= 1) & (n2 <= ntemp)
        cbind(n1[good], n2[good])
    }
    if(zero.out) dat[dat==0] <- NA
    N <- ncol(dat)
    M <- nrow(dat)
    m <- min(c(round(R/dx), M))
    n <- min(c(round(R/dy), N))
    ind <- rbind(as.matrix(expand.grid(0, 1:n)), as.matrix(expand.grid(1:m,
        0)), as.matrix(expand.grid(c(-(m:1), 1:m), 1:n)))
    d <- sqrt((dx * ind[, 1])^2 + (dy * ind[, 2])^2)
    good <- (d > 0) & (d <= R)
    ind <- ind[good, ]
    d <- d[good]
    ind <- ind[order(d), ]
    d <- sort(d)
    nbin <- nrow(ind)
    holdVG <- rep(NA, nbin)
    holdN <- rep(NA, nbin)
    for (k in 1:nbin) {
        MM <- SI(M, ind[k, 1])
        NN <- SI(N, ind[k, 2])
        numNA <- sum(is.na(dat[MM[,1],NN[,1]]) | is.na(dat[MM[,2],NN[,2]]),na.rm=TRUE)
        holdN[k] <- length(MM) * length(NN) - numNA
        BigDiff <- (dat[MM[, 1], NN[, 1]] - dat[MM[, 2], NN[,2]])
        holdVG[k] <- mean(0.5 * (BigDiff)^q, na.rm=TRUE)
    }
    top <- tapply(holdVG * holdN, d, FUN = "sum")
    bottom <- tapply(holdN, d, FUN = "sum")
    dcollapsed <- as.numeric(names(bottom))
    vgram <- top/bottom
    dimnames(vgram) <- NULL
    out <- list(vgram = vgram, d = dcollapsed, ind = ind, d.full = d,
        vgram.full = holdVG, N = holdN, dx = dx, dy = dy, q=q)
    class(out) <- "structurogram.matrix"
    return(out)
} # end of 'structurogram.matrix' function.

plot.structurogram.matrix <- function(x,...) {
   par(mfrow=c(1,2), bg="beige")
   plot( x$d, x$vgram, xlab="separation distance", ylab=paste("structure (q=", x$q, ")", sep=""), ...)
   points( x$d.full, x$vgram.full, pch=".")
   plot.vgram.matrix(x, main="Structure by direction")
   invisible()
} # end of 'plot.structurogram.matrix' function.

structurogram <- function(loc, y, q=2, id = NULL, d = NULL, lon.lat = FALSE, dmax = NULL, N = NULL, breaks = NULL) 
{
    y <- cbind(y)
    if (is.null(id)) {
        n <- nrow(loc)
        ind <- rep(1:n, n) > rep(1:n, rep(n, n))
        id <- cbind(rep(1:n, n), rep(1:n, rep(n, n)))[ind, ]
    }
    if (is.null(d)) {
        loc <- as.matrix(loc)
        if (lon.lat) {
            d <- rdist.earth(loc)[id]
        }
        else {
            d <- rdist(loc, loc)[id]
        }
    }
    vg <- 0.5 * rowMeans(cbind((y[id[, 1], ] - y[id[, 2], ])^q), 
        na.rm = TRUE)
    call <- match.call()
    if (is.null(dmax)) {
        dmax <- max(d)
    }
    od <- order(d)
    d <- d[od]
    vg <- vg[od]
    ind <- d <= dmax & !is.na(vg)
    out <- list(d = d[ind], val = vg[ind], call = call, q=q)
    if (!is.null(breaks) | !is.null(N)) {
        out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks))
    }
    class(out) <- "structurogram"
    return(out)
} # end of 'structurogram' function.

plot.structurogram <- function(x,...) {
   plot(x$d, x$val, xlab="separation distance", ylab=paste("structure (q=",x$q,")", sep=""), ...)
   lines(x$centers, x$stats["mean",], col="darkblue")
   invisible()
} # end of 'plot.structurogram' function.
