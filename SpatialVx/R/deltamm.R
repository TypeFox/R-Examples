deltamm <- function(x, p = 2, max.delta = Inf, const = Inf, verbose = FALSE, ...) {

    if(verbose) begin.tiid <- Sys.time()

    out <- x
    a <- attributes(x)

    out$match.type <- "deltamm"
    out$match.message <- "Objects merged/matched based on the Baddeley Delta metric via the deltamm function."

    # below is an internal function to calculate Baddeley's delta metric for each pair
    # of features (possibly already merged) as indicated by 'id', which gives the observed
    # feature number in its first position and the forecast feature number in the second.

    bfun <- function(id, OB, FC, p, const, verbose) {

	# 'OB' is a list of the observed distance maps.
	# 'FC' is a list of the forecast distance maps.

	# 'p', 'const', and 'verbose' are same as arguments passed to 'deltamm'.

        j <- id[1]
        k <- id[2]

	if(verbose) cat("Calculating Baddeley Delta Metric between forecast feature ", k, " and observed feature ", j, "\n")

	# First, try to figure out if one of 'OB' or 'FC' is a list of lists or just a list.

	if(class(OB[[ 1 ]]) == "im") {

	    A <- OB[[ j ]]
            if(class(FC[[ 1 ]]) == "im") B <- FC[[ k ]]
	    else B <- FC[[ j ]][[ k ]]

	} else {

	    A <- OB[[ k ]][[ j ]]
	    B <- FC[[ k ]]

	}

        if(!is.infinite(const)) {

	    A <- eval.im(pmin.int(A, const))
	    B <- eval.im(pmin.int(B, const))

	} # end of if 'const' is infinite stmt.

	if(is.infinite(p)) {

	    res <- eval.im(abs(A - B))
	    res <- summary(res)$max

	} else {

	    res <- eval.im(abs(A - B)^p)
	    res <- summary(res)$mean
	    res <- res^(1/p)

	} # end of if else 'p' is infinite stmts.

	return(res)

    } # end of internal 'bfun' function.

    # If this is the first pass, do the first pass.  If it is a subsequent pass, must first figure out what
    # has already been matched/merged, and then proceed.


    X <- x$X.feats
    Xhat <- x$Y.feats 

    xdim <- dim(X[[1]][["m"]])

    n <- length(X)
    m <- length(Xhat)

    no.X <- is.null(X) || n == 0
    no.Xhat <- is.null(Xhat) || m == 0

    if(no.X || no.Xhat) {

        # No features in or both fields.  Return an appropriate object.

	if(no.X && no.Xhat && verbose) cat("No identified features in either field.  Therefore, no matches/merges made.\n")
	else if(no.X && verbose) cat("No identified observed features.  Therefore, no matches/merges made.\n")
	else if(no.Xhat && verbose) cat("No identified model features.  Therefore, no matches/merges made.\n")

	funmatched <- 1:m
	vxunmatched <- 1:n
	matches <- cbind(integer(0), integer(0))
	merges <- NULL

	out$unmatched <- list(X = vxunmatched, Xhat = funmatched)
	out$matches <- matches
  	out$merges <- merges

	class(out) <- "matched"
	return(out)

    } # end of if no features in either or both field stmts.

    # In order to reduce the number of computations, first find distance maps for
    # all individual features so that they can be re-used when necessary.  To do so,
    # need to find a bounding box for the entire set of features.

    # Find the bounding box for the union of all of the feature fields.
    # For some reason, this is not working with boundingbox directly.
    # So, doing it manually.

    if(verbose) cat("Finding a bounding box for the entire set of features.\n")

    Xbb <- do.call( boundingbox, X )
    Ybb <- do.call( boundingbox, Xhat )

    bb <- boundingbox(Xbb, Ybb)

    if(verbose) {

        cat("Finished finding overall bounding box for the entire set of features.\n")
        print(bb)

    } # end of if 'verbose' stmt.

    Xbox <- lapply(X, rebound, rect = bb)
    Ybox <- lapply(Xhat, rebound, rect = bb)

    if(verbose) cat("Finding the distance maps for each feature in each field.\n")

    dX <- lapply(Xbox, distmap, ...)
    dXhat <- lapply(Ybox, distmap, ...)

    # Ready for step 1.
    if(verbose) cat("Step 1: Finding Upsilon matrix containing the Baddeley delta between each individual feature across fields.\n")

    ind <- cbind(rep(1:n, m), rep(1:m, each = n))

    Upsilon <- apply(ind, 1, bfun, OB = dX, FC = dXhat, p = p, const = const, verbose = verbose)

    Upsilon <- matrix(Upsilon, n, m)
    # Upsilon should be a matrix of delta metrics whose rows
    # are observed feature numbers and whose columns are
    # forecast feature numbers so that the (j, k)-th element
    # is the metric between observed feature j and forecast feature k.

    if(verbose) {

        cat("Step 1 completed.  Upsilon matrix given by:\n")
        print(Upsilon)

    } # end of if 'verbose' stmt.

    Psi <- Ksi <- matrix(NA, n, m)
    # Ksi is j-th observed compared with forecast mergings
    # Psi is k-th forecast compared with observed mergings
    
    # Note that I've transposed the matrices from those in the paper.

    # Determine the rank order of each row and column of
    # Upsilon.  Ensure result is an n by m matrix.
    o.Ksi <- t(apply(Upsilon, 1, order, na.last = TRUE))
    o.Psi <- apply(Upsilon, 2, order, na.last = TRUE)

    # Take the first column of Ksi and first row of Psi
    # directly from Upsilon.
    Ksi[,1] <- apply(Upsilon, 1, min, na.rm = TRUE)
    Psi[1,] <- apply(Upsilon, 2, min, na.rm = TRUE)

    # Now compute the distance maps for merged features.
    if(verbose) cat("Finding distance maps for potential merges.\n")

    if(verbose) cat(" Finding ", m - 1, " potential forecast merges for each of the ", n, " observed features.\n")

    # Ksi.dmap and Psi.dmap are lists of lists
    # that hold the distance maps for merges.

    Ksi.dmap <- Psi.dmap <- list()

    for(j in 1:n) {

        if(verbose) cat(j, " ")

        o <- (1:m)[ o.Ksi[j,] ]

        # The first "union" is just the individual feature with
        # minimum delta, so start loop at 2.  Also, the last
        # merge is always the union of all features, so need only
        # compute that distance map once.  To be done later...

        newobj <- list()

        for(i in 2:(m - 1)) newobj[[ i - 1 ]] <- union.owin(Xhat[[ o[ i - 1] ]], Xhat[[ o[ i ] ]])

        newobj2 <- lapply(newobj, rebound, rect = bb)

        Ksi.dmap[[ j ]] <- lapply(newobj2, distmap, ...)

        # Find the distance map for the merging of all features.
        # This is the same for all n cases, so need only be done once.

        if(j == n) {

	    newlast <- union.owin(newobj[[ m - 2 ]], Xhat[[ o[ m ] ]])
	    newlast <- rebound(newlast, bb)
	    Ksi.last <- distmap(newlast, ...)

        } # end of if last iteration of loop stmt.

    } # end of outer 'j' loop.


    if(verbose) cat("\nFinding ", n - 1, " potential observed merges for each of the ", m, " forecast features.\n")

    for(k in 1:m) {

        if(verbose) cat(k, " ")

        o <- (1:n)[ o.Psi[, k] ]

        newobj <- list()

        for(i in 2:(n - 1)) newobj[[ i - 1 ]] <- union.owin(X[[ o[ i - 1] ]], X[[ o[ i ] ]])

        newobj2 <- lapply(newobj, rebound, rect = bb)

        Psi.dmap[[ k ]] <- lapply(newobj2, distmap, ...)

        if(k == m) {

	    newlast <- union.owin(newobj[[ n - 2 ]], X[[ o[ n ] ]])
	    newlast <- rebound(newlast, bb)
	    Psi.last <- distmap(newlast, ...)

        }

    } # end of outer 'k' loop.

    if(verbose) cat("\nAll necessary distance maps have been found.  Calculating delta metrics for forecast merges.\n")

    # Grab the indicators for the next round but leave out metrics that have already been calculated.
    ind1 <- ind2 <- ind

    ind1[ ind[,2] == 1,] <- NA
    ind1[ ind[,2] == m,] <- NA
    ind1 <- ind1[!is.na(ind1[,1]),]
    ind1[,2] <- ind1[,2] - 1

    ind2[ ind[,1] == 1,] <- NA
    ind2[ ind[,1] == n,] <- NA
    ind2 <- ind2[!is.na(ind2[,1]),]
    ind2[,1] <- ind2[,1] - 1

    tmp <- apply(ind1, 1, bfun, OB = dX, FC = Ksi.dmap, p = p, const = const, verbose = verbose)

    Ksi[, 2:(m - 1)] <- tmp

    ind1 <- cbind(1:n, rep(1, n))
    Ksi.last2 <- list()
    Ksi.last2[[ 1 ]] <- Ksi.last

    tmp <- apply(ind1, 1, bfun, OB = dX, FC = Ksi.last2, p = p, const = const, verbose = verbose)

    Ksi[, m] <- tmp

    if(verbose) cat("\nAll metrics for forecast merges found.  Calculating delta metrics for observed merges.\n")

    tmp <- apply(ind2, 1, bfun, OB = Psi.dmap, FC = dXhat, p = p, const = const, verbose = verbose)

    Psi <- t(Psi)
    Psi[,2:(n - 1)] <- tmp

    ind2 <- cbind(rep(1, m), 1:m)
    Psi.last2 <- list()
    Psi.last2[[ 1 ]] <- Psi.last

    tmp <- apply(ind2, 1, bfun, OB = Psi.last2, FC = dXhat, p = p, const = const, verbose = verbose)

    Psi[, n] <- tmp
    Psi <- t(Psi)

    bigQ <- array(NA, dim = c(n, m, 3))

    bigQ[,,1] <- Upsilon
    bigQ[,,2] <- Psi
    bigQ[,,3] <- Ksi

    out$Q <- bigQ

    if(all(bigQ > max.delta)) {

	# Case of no matches or merges.  Return an appropriate object.

	if(verbose) cat("\nAll delta metrics are larger than max.delta, no merges/matches.\n")

	funmatched <- 1:m
        vxunmatched <- 1:n
        matches <- cbind(integer(0), integer(0))
        merges <- NULL

        out$unmatched <- list(X = vxunmatched, Xhat = funmatched)
        out$matches <- matches
	out$merges <- merges

        class(out) <- "matched"
        return(out)

    } else {

	if(verbose) {

	    cat("Psi:\n")
	    print(Psi)

	    cat("Ksi:\n")
	    print(Ksi)

	    cat("\nAll Baddeley metrics found.  Book keeping ...\n")

	} # end of if 'verbose' stmt.

        # J is a list of integer vectors with observed single and merged features.
        # K is the same but for forecast features.
        J <- K <- list()

        ind <- cbind(rep(1:n, 3 * m), rep(rep(1:m, each = n), 3), rep(1:3, each = n * m), 1:(3 * n * m))

	for(jk in 1:(3 * n * m)) {

	   if(jk <= n * m) {

		# Case of Upsilon matrix.

		J[[ jk ]] <- ind[ jk, 1 ]
		K[[ jk ]] <- ind[ jk, 2 ]

	    } else if((jk > n * m) && (jk <= 2 * n * m)) {

		# Case of 'Psi' matrix (observed mergings).

		K[[ jk ]] <- ind[ jk, 2 ]

		jj <- 1:(ind[ jk, 1])
		J[[ jk ]] <- (1:n)[ o.Psi[ jj, ind[jk, 2] ] ]

	    } else {

		# Case of 'Ksi' matrix (forecast mergings).

		J[[ jk ]] <- ind[ jk, 1 ]

		kk <- 1:(ind[jk, 2])
		K[[ jk ]] <- (1:m)[ o.Ksi[ ind[jk, 1], kk ] ]

	    } # end of if else which matrix are we on stmts.

	} # end of for 'jk' loop.

	iter <- 1
	matches <- cbind(integer(0), integer(0))

	nn <- 1:n
        mm <- 1:m

	bigQ[ bigQ > max.delta ] <- NA
	bigQ <- c(bigQ)

	efun <- function(x, ftr) return(any(is.element(ftr, x)))

	while( length(nn) > 0 && length(mm) > 0 && iter < 3 * n * m + 1) {

	    if(any(!is.na(bigQ))) {

	        minQ <- min(bigQ, na.rm = TRUE)
	        id <- (ind[,4])[ bigQ == minQ ]
		id <- id[ !is.na(id) ]
		id <- id[ 1 ]

		vx <- J[[ id ]]
		fc <- K[[ id ]]

		newmatch <- cbind(fc, vx)
		matches <- rbind(matches, newmatch)

		# 'matches' has been updated, now need to eliminate
		# features that have been matched.

		bigQ[id] <- NA
		id1 <- unlist(lapply(J, efun, ftr = vx))
		id2 <- unlist(lapply(K, efun, ftr = fc))

		bigQ[ id1 | id2 ] <- NA
		nn <- nn[ !is.element(nn, vx) ]
		mm <- mm[ !is.element(mm, fc) ]
		iter <- iter + 1

	    } else iter <- Inf

	} # end of while there are still possible matches/merges loop.

	out$unmatched <- list(X = nn, Xhat = mm)

	colnames(matches) <- c("Forecast", "Observed")
	out$matches <- matches

	merges <- MergeIdentifier(matches)
	out$merges <- merges

    } # end of if no deltas small enough to match/merge stmts.

    if(verbose) print(Sys.time() - begin.tiid)

    class(out) <- "matched"
    return(out)

} # end of 'deltamm' function.
