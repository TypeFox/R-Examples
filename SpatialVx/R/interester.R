interester <- function(x, properties = c("cent.dist", "angle.diff", "area.ratio", 
    "int.area", "bdelta", "haus", "ph", "med", "msd", 
    "fom", "minsep"), weights = c(0.24, 0.12, 0.17, 0.12, 0, 0, 0, 0, 0, 0, 0.35),
    b1 = c(35, 30, 0, 0, 0.5, 35, 20, 40, 120, 1, 40),
    b2 = c(100, 90, 0.8, 0.25, 85, 400, 200, 200, 400, 0.25, 200), 
    verbose = FALSE, ...) {

    if(verbose) begin.tiid <- Sys.time()

    if((is.null(x$X.feats) || length(x$X.feats) == 0) || (is.null(x$Y.feats) || length(x$Y.feats) == 0)) {

	warning("interester: No features in one or both of the fields.  Returning NULL.")
	return(NULL)

    }

    zerow <- weights == 0
    np <- sum(!zerow)

    if(all(zerow)) {

	warning("interester: all weights are zero so that no interest is calculated.  Returning NULL.")
	return(NULL)

    }

    if(any(zerow)) {

	properties <- properties[ !zerow ]
	weights <- weights[ !zerow ]
	b1 <- b1[ !zerow ]
	b2 <- b2[ !zerow ]

    }

    # Find slope and intercept terms based on property type, b0 and b1.
    a0 <- a1 <- numeric(np) + NA

    type.ind <- !is.element(properties, c("area.ratio", "int.area", "fom"))

    if(any(type.ind)) {

	a1[ type.ind ] <- -1/(b2[ type.ind ] - b1[ type.ind ])
	a0[ type.ind ] <- 1 - a1[ type.ind ] * b1[ type.ind ]

    }

    type.ind <- is.element(properties, c("area.ratio", "int.area"))

    if(any(type.ind)) {

	a1[ type.ind ] <- 1/(b2[ type.ind ] - b1[ type.ind ])
	a0[ type.ind ] <- 1 - a1[ type.ind ] * b2[ type.ind ]

    }

    ipwlin <- function(x, b1, b2, a0, a1, property, ...) {

	dn <- is.element(property, c("cent.dist", "angle.diff", "bdelta", "haus", "ph", "med", "msd", "minsep"))
	up <- is.element(property, c("area.ratio", "int.area"))
	fom <- property == "fom"

	res <- numeric(length(x))

	if(any(is.na(x))) {

	    if( any(dn & is.na(x)) ) x[ dn & is.na(x) ] <- b2[ dn & is.na(x) ] + 1e-8

	    if( any(up & is.na(x)) ) x[ up & is.na(x) ] <- b1[ up & is.na(x) ] - 1e-8

	    if( any(fom) ) x[ fom & is.na(x) ] <- 100

	} # end of if any missing x values stmts.

        if(any(dn)) {

	    dn2 <- x <= b1
	    dn3 <- x > b2

	    if(any(dn & dn2)) res[ dn & dn2 ] <- 1
	    if(any(dn & dn3)) res[ dn & dn3 ] <- 0
	    if(any(dn & !(dn2 | dn3))) {

		look <- a0[ dn & !(dn2 | dn3) ] + a1[ dn & !(dn2 | dn3) ] * x[ dn & !(dn2 | dn3) ]
		if(any(look < 0)) look[ look < 0 ] <- 0
		res[ dn & !(dn2 | dn3) ] <- look

	    }

        } 

 	if(any(up)) {

	    up2 <- x < b1
	    up3 <- x >= b2

	    if(any(up & up2)) res[ up & up2 ] <- 0
	    if(any(up & up3)) res[ up & up3 ] <- 1
	    if(any(up & !(up2 | up3))) {

		look <- a0[ up & !(up2 | up3) ] + a1[ up & !(up2 | up3) ] * x[ up & !(up2 | up3) ]
		if(any(look > 1)) look[ look > 1 ] <- 1
		res[ up & !(up2 | up3) ] <- look

	    }

        } 

	if(any(fom)) {

	    look <- b1[ fom ] * exp( -0.5 * ((x[ fom ] - 1) / b2[ fom ])^4 )
	    if(look < 0) look <- 0
	    else if(look > 1) look <- 1
            res[ fom ] <- look

        }

	names(res) <- names(x)
        return(res)

    } # end of internal 'ipwlin' function.

    ifun <- function(id, bigX, Xhat, b1, b2, a0, a1, p, ...) {

	if(verbose) cat("Forecast feature: ", id[2], " vs Observed feature: ", id[1], "\n")

	Xtmp <- bigX[[ id[1] ]]
	Ytmp <- Xhat[[ id[2] ]]

	A <- FeatureComps(Ytmp, Xtmp, which.comps = p, ...)

	A <- distill(A)

	nomen <- names(A)

	res <- ipwlin(A, b1 = b1, b2 = b2, a0 = a0, a1 = a1, property = p, ...)

	if(is.element("angle.diff", p)) {

	    aspX <- FeatureAxis(Xtmp)$aspect.ratio
	    aspY <- FeatureAxis(Ytmp)$aspect.ratio
	    conX <- ((aspX - 1)^2 / (aspX^2 + 1))^(0.3)
	    conY <- ((aspY - 1)^2 / (aspY^2 + 1))^(0.3)
	    con <- sqrt( conX * conY )
	    if(length(con) == 0) con <- 0
	    res[ "angle.diff" ] <- res[ "angle.diff" ] * con

	}

	if(all(is.element(c("cent.dist", "area.ratio"), p)) ) {

	    res[ "cent.dist" ] <- res[ "cent.dist" ] * res[ "area.ratio" ]

	}

	return(res)

    } # end of internal 'ifun' function.

    N <- length(x$X.feats)
    M <- length(x$Y.feats)
    ind <- cbind(rep(1:N, M), rep(1:M, each = N))

    if(verbose) cat("\n\nFinding interest between each pair of ", N, " observed features and ", M, " forecast features.\n\n")
    res1 <- apply(ind, 1, ifun, bigX = x$X.feats, Xhat = x$Y.feats, b1 = b1, b2 = b2, a0 = a0, a1 = a1, p = properties, ...)

    out <- list()
    a <- attributes(x)
    a$names <- NULL
    attributes(out) <- a

    out$interest <- res1
    out$total.interest <- matrix( colSums(res1 * matrix(weights, np, N * M), na.rm = TRUE), N, M)

    if(verbose) print(Sys.time() - begin.tiid)

    class(out) <- "interester"
    return(out)

} # end of 'interester' function.

print.interester <- function(x, ...) {

    i <- x$interest
    if(prod(dim(i)) < 100) {

	cat("Pair-by-pair interest values:\n")
	print(i)

    } else {

	cat("Summary of pair-by-pair interest values.\n")
	print(summary(t(i)))

    }

    cat("Total interest (pair-by-pair)\n")
    print(x$total.interest)

    invisible()

} # end of 'print.interester' function.

summary.interester <- function(object, ..., min.interest = 0.8, long = TRUE, silent = FALSE) {

    x <- object
    a <- attributes(x)
    a$names <- NULL

    out <- list()
    attributes(out) <- a

    i <- x$interest
    ti <- x$total.interest
    d <- dim(ti)

    ind <- cbind(rep(1:d[1], d[2]), rep(1:d[2], each = d[1]))

    y <- cbind(ind, c(ti))
    o <- order(c(ti), decreasing = TRUE)
    y <- y[o,]
    colnames(y) <- c("obs feature", "mod feature", "total interest")

    if(!silent) {

	cat("Ranking of feature pairings based on total interest.\n")
	id <- y[,3] >= min.interest
	if(any(id)) print(y[id,])
	if(any(id) && any(!id) && long) cat("\n--------\n")
	if(any(!id) && long) print(y[!id,])

    }

    out$sorted.interest <- y

    if(d[1] > 1 && d[2] > 1) mmi <- quantile( c(c(apply(ti, 1, max, na.rm = TRUE)), c(apply(ti, 2, max, na.rm = TRUE))), probs = 0.5)
    else if((d[1] > 1 && d[2] == 1) || (d[1] == 1 && d[2] > 1) || (is.null(d) && length(ti) > 0)) mmi <- max(c(ti), na.rm = TRUE)
    else mmi <- NULL

    if(!silent) {

	if(d[1] > 1 && d[2] > 1) cat("Median of Maximum Interest = ", mmi, "\n")
	else if((d[1] > 1 && d[2] == 1) || (d[1] == 1 && d[2] > 1) || (is.null(d) && length(ti) > 0)) cat("Maximum Interest = ", mmi, "\n")

    }

    if(d[1] > 1 && d[2] > 1) mmi.msg <- paste("Median of Maximum Interest = ", mmi, sep = "")
    else if((d[1] > 1 && d[2] == 1) || (d[1] == 1 && d[2] > 1) || (is.null(d) && length(ti) > 0)) mmi.msg <- paste("Maximum Interest = ", mmi, sep = "")

    out$message <- mmi.msg
    out$mmi <- mmi 

    class(out) <- "summary.interester"

    invisible(out)

} # end of 'summary.interester' function.

print.summary.interester <- function(x, ..., min.interest = 0.8, long = TRUE) {

    y <- x$sorted.interest
    cat("Ranking of feature pairings based on total interest.\n")
    id <- y[,3] >= min.interest
    if(any(id)) print(y[id,])
    if(any(id) && any(!id) && long) cat("\n--------\n")
    if(any(!id) && long) print(y[!id,])

    cat(x$message, "\n")
    print(x$mmi)

    invisible()

} # end of 'print.summary.interester' function.

distill.FeatureComps <- function(x, ...) {

    desired.names <- c("cent.dist", "angle.diff", "area.ratio", 
        "int.area", "bdelta", "haus", "ph", "med", "msd", 
        "fom", "minsep")

    nomen <- names(x)

    id <- is.element(nomen, desired.names)

    # Grab the relevant names and re-order them.
    nomen <- nomen[ id ]
    nomen <- desired.names[ is.element(desired.names, nomen) ]

    out <- numeric(0)
    nomenout <- character(0) 

    for( i in 1:length(nomen)) {

	tmp <- x[[ nomen[ i ] ]]

	tdim <- dim(tmp)

	if(!is.null(tdim)) {

	    n <- length(c(tmp))
	    if(n == 0) tmp <- NA
	    out <- c(out, c(tmp))
	    tnames <- paste(rep(1:tdim[1], tdim[2]), rep(1:tdim[2], each = tdim[1]), sep = "")

	    if(n > 1) nomenout <- c(nomenout, paste(nomen[ i ], tnames, sep = ""))
	    else nomenout <- c(nomenout, nomen[ i ])

	} else {

	    n <- length(tmp)
	    if(n == 0) tmp <- NA
	    out <- c(out, tmp)
	    if(n > 1) nomenout <- c(nomenout, paste(nomen[ i ], 1:n, sep = ""))
	    else nomenout <- c(nomenout, nomen[ i ])

	}

    } # end of for 'i' loop.

    names(out) <- nomenout

    return(out)

} # end of 'distill.FeatureComps' function.
