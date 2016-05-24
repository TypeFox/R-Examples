FeatureFinder <- function (object, smoothfun = "disk2dsmooth", do.smooth = TRUE, smoothpar = 1, 
    smoothfunargs = NULL, thresh = 1e-08, idfun = "disjointer", min.size = 1, max.size = Inf,
    fac = 1, zero.down = FALSE, time.point = 1, model = 1, ...) 

{

    theCall <- match.call()

    if( length( min.size ) == 1 ) min.size <- rep( min.size, 2 )
    if( length( max.size ) == 1 ) max.size <- rep( max.size, 2 )

    if( length( min.size ) != 2 ) stop("FeatureFinder: invalid min.size argument.  Must have length one or two.")
    if( length( max.size ) != 2 ) stop("FeatureFinder: invalid max.size argument.  Must have length one or two.")

    if( any( min.size < 1 ) ) stop("FeatureFinder: invalid min.size argument.  Must be >= 1.")
    if( any( max.size < min.size ) ) stop("FeatureFinder: invalid max.size argument.  Must be >= min.size argument.")

    a <- attributes(object)

    if ( !missing(time.point) && !missing(model) ) dat <- datagrabber( object, time.point = time.point, model = model )
    else if (!missing(time.point)) dat <- datagrabber( object, time.point = time.point )
    else if (!missing(model)) dat <- datagrabber( object, model = model )
    else dat <- datagrabber( object )

    X <- dat$X
    Y <- dat$Xhat

    xdim <- a$xdim

    if( do.smooth ) {

	if( length( smoothpar ) == 1 ) smoothpar <- rep( smoothpar, 2 )
	else if( length( smoothpar ) > 2 ) stop("FeatureFinder: invalid smoothpar argument.  Must have length one or two.")

        Xsm <- do.call(smoothfun, c(list(x = X, lambda = smoothpar[ 2 ]), smoothfunargs))
        Ysm <- do.call(smoothfun, c(list(x = Y, lambda = smoothpar[ 1 ]), smoothfunargs))

	if (zero.down) {

            Xsm[ Xsm < 0 ] <- 0
            Xsm <- zapsmall( Xsm )
            Ysm[ Ysm < 0 ] <- 0
            Ysm <- zapsmall( Ysm )

        }

    } else {

	Xsm <- X
	Ysm <- Y

    } # end of if 'do.smooth' stmts.

    # if (is.null(thresh)) thresh <- c(quantile(c(Xsm), 0.75), quantile(c(Ysm), 0.75))
    # else if (length(thresh) == 1) thresh <- c(thresh, thresh)

    if (length(thresh) == 1) thresh <- c(thresh, thresh)

    thresh <- thresh * fac

    sIx <- sIy <- matrix(0, xdim[1], xdim[2])

    sIx[ Xsm >= thresh[1] ] <- 1
    sIy[ Ysm >= thresh[2] ] <- 1

    X.feats <- do.call(idfun, c(list(x = sIx), list(...)))
    Y.feats <- do.call(idfun, c(list(x = sIy), list(...)))

    # Check if any features exist.
    if( length( X.feats ) == 0 ) X.feats <- NULL
    if( length( Y.feats ) == 0 ) Y.feats <- NULL

    # Now remove any features that are too small or possibly too big.
    if( any( min.size > 1 ) | any( max.size < prod( xdim ) ) ) {

	Nfun <- function(Obj) return( sum( colSums( as.matrix( Obj ), na.rm = TRUE ), na.rm = TRUE ) )

	if( !is.null( X.feats ) ) {

            Xnums <- c(unlist(lapply(X.feats, Nfun)))
	    Xnums <- Xnums >= min.size[ 2 ] & Xnums <= max.size[ 2 ]
	    Xj <- ( 1:length( Xnums ) )[ Xnums ]

	    X.feats0 <- list()

	    if( length( Xj ) > 0 ) {

		for( i in 1:length( Xj ) ) X.feats0[[ i ]] <- X.feats[[ Xj[ i ] ]]
		X.feats <- X.feats0

	    }

	} # end of if X features are present stmt.

	if( !is.null( Y.feats ) ) {

            Ynums <- c( unlist( lapply( Y.feats, Nfun ) ) )
            Ynums <- Ynums >= min.size[ 1 ] & Ynums <= max.size[ 1 ] 
            Yj <- (1:length(Ynums))[ Ynums ]

	    Y.feats0 <- list()

	    if( length( Yj ) > 0 ) {

		for( i in 1:length( Yj ) ) Y.feats0[[ i ]] <- Y.feats[[ Yj[ i ] ]]
		Y.feats <- Y.feats0

	    }

	} # end of if Y features are present stmt.

    } # end of if remove features that are too big or two small.

    Xlab <- Ylab <- matrix(0, xdim[1], xdim[2])

    if (!is.null(X.feats)) for (i in 1:length(X.feats)) Xlab[ as.matrix( X.feats[[i]] ) ] <- i
    else X.feats <- NULL

    if (!is.null(Y.feats)) for (j in 1:length(Y.feats)) Ylab[ as.matrix( Y.feats[[j]] ) ] <- j
    else Y.feats <- NULL

    out <- list()

    attributes(out) <- a

    out$X <- X
    out$Xhat <- Y
    out$X.feats <- X.feats
    out$Y.feats <- Y.feats
    out$X.labeled <- Xlab
    out$Y.labeled <- Ylab
    out$identifier.function <- "convthresh"
    out$identifier.label <- "Convolution Threshold"

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model

    if (length(a$data.name) == a$nforecast + 2) {

        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[1:2]

    } else {

        dn <- a$data.name[-1]
        vxname <- a$data.name[1]

    }

    if (!is.numeric(model)) 
        model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    attr(out, "data.name") <- c(vxname, dn[model.num])
    attr(out, "call") <- theCall
    class(out) <- "features"

    return(out)

} # end of 'FeatureFinder' function.
