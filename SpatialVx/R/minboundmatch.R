minboundmatch <- function(x, type = c("single", "multiple"), mindist = Inf, verbose = FALSE, ...) {

    if(verbose) begin.tiid <- Sys.time()

    if(class( x ) != "features") stop("minboundmatch: invalid x argument.")

    type <- tolower( type )
    type <- match.arg( type )

    out <- x
    a <- attributes( x )

    out$match.type <- "minboundmatch"
    out$match.message <- paste("Matching based on minimum boundary separation using ", type, " matches.", sep = "")

    Xfeats <- x$X.feats
    Yfeats <- x$Y.feats

    if( !is.null(Xfeats) ) n <- length( Xfeats )
    else n <- 0

    if( !is.null(Yfeats) ) m <- length( Yfeats )
    else m <- 0

    if(m == 0 && n == 0) {

	if(verbose) cat("\n", "No features detected in either field.  Returning NULL.\n")
	return(NULL)

    } else if(m == 0) {

	if(verbose) cat("\n", "No features detected in forecast field.  Returning NULL.\n")
        return(NULL)

    } else if(n == 0) {

	if(verbose) cat("\n", "No features detected in observed field.  Returning NULL.\n")
        return(NULL)

    } # end of quietly return NULL if no features in one or both fields stmts.

    ind <- cbind( rep(1:n, m), rep(1:m, each = n) )

    Xdmaps <- lapply(Xfeats, distmap, ...)

    minsepfun <- function(id, dm, Y) {

	i <- id[ 1 ]
	j <- id[ 2 ]

	X <- dm[[ i ]]
	X <- X$v
	Ym <- as.matrix( Y[[ j ]] )

	Z <- X[ Ym ]

	return( min( c( Z ), na.rm = TRUE) )

    } # end of 'minsep' internal function.

    res <- apply(ind, 1, minsepfun, dm = Xdmaps, Y = Yfeats)

    res <- cbind(ind, res)
    colnames( res ) <- c( "Observed Feature No.", "Forecast Feature No.", "Minimum Boundary Separation" )
    good <- res[, 3] <= mindist
    res <- res[ good, , drop = FALSE]

    out$values <- res

    # Above is a simple numeric vector of minimum boundary separation numbers.

    o <- order( res[, 3] )
    res <- res[o, , drop = FALSE]

    if(type == "single") {

	N <- dim( res )[ 1 ]
	id <- 1:N

	id <- id[ o ]

	matches <- cbind( numeric( 0 ), numeric( 0 ) )

	for(i in 1:N) {

	    matches <- rbind( matches, res[1, 2:1 ] )
	    id2 <- (res[, 1] == res[1, 1]) | (res[, 2] == res[1, 2])

	    res <- res[ !id2, , drop = FALSE ]
	    id <- id[ !id2 ]

	    if(length( id ) == 0) break

	} # end of for 'i' loop.

    } else {

	matches <- res[, 2:1 , drop = FALSE]

	matchlen <- dim( matches )[ 1 ]
        fuq <- unique( matches[, 1 ] )
        flen <- length( fuq )

        ouq <- unique( matches[, 2 ] )
        olen <- length( ouq )

	if(matchlen > 0) {

            if(matchlen == flen && matchlen > olen) {

                if(verbose) cat("Multiple observed features are matched to one or more forecast feature(s).  Determining implicit merges.\n")

            } else if(matchlen > flen && matchlen == olen) {

                if(verbose) cat("Multiple forecast features are matched to one or more observed feature(s).  Determining implicit merges.\n")

            } else if(matchlen > flen && matchlen > olen) {

                if(verbose) cat("Multiple matches have been found between features in each field.  Determining implicit merges.\n")

            } else if(matchlen == flen && matchlen == olen) {

                if(verbose) cat("No multiple matches were found.  Thus, no implicit merges need be considered.\n")

            } # end of if else which fields have multiple matches stmts.

            out$implicit.merges <- MergeIdentifier( matches )

            } else {

                if(verbose) cat("No objects matched.\n")
                out$implicit.merges <- NULL

            } # end of if any matches stmts.

    } # end of if else type is single or multiple stmts.

    matches <- matches[ order( matches[, 1] ), , drop = FALSE]

    colnames( matches ) <- c( "Forecast", "Observed" )
    out$matches <- matches

    out$unmatched <- list( X = unique( ind[, 1] )[ !is.element( unique(ind[, 1]), matches[, 2] ) ],
                                Xhat = unique( ind[, 2] )[ !is.element( unique( ind[, 2] ), matches[, 1] ) ] )

    if(verbose) print(Sys.time() - begin.tiid)

    class( out ) <- "matched"
    return( out )

} # end of 'minboundmatch' function.
