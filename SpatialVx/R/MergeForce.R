MergeForce <- function(x, verbose = FALSE) {

    if(class(x) != "matched") stop("MergeForce: invalid x type.")

    a <- attributes(x)
    out <- list()
    a$names <- NULL
    attributes(out) <- a

    if(!is.null(x$implicit.merges) && !is.null(x$merges)) {

	warning("MergeForce: both implicit.merges and merges components found.  Using merges component.")
	m <- x$merges

    } else if(is.null(x$implicit.merges) && is.null(x$merges)) {

	if(verbose) cat("\nNo merges found.  Returning object x unaltered.\n")
	return(x)

    } else if(is.null(x$merges)) {

	m <- x$implicit.merges

    } else {

	m <- x$merges

    }

    out$X <- x$X
    out$Xhat <- x$Xhat
    out$identifier.function <- x$identifier.function
    out$identifier.label <- x$identifier.label
    out$match.type <- c(x$match.type, "MergeForce") # TO DO: make sure this does not cause problems.
    out$match.message <- paste(x$match.message, " (merged) ", sep = "")

    maxF <- max(x$Y.labeled, na.rm = TRUE)
    maxO <- max(x$X.labeled, na.rm = TRUE)

    xdim <- dim(x$X.labeled)

    # Number of matches (after merging)
    nmatches <- length(m)

    matches <- cbind(1:nmatches, 1:nmatches)
    colnames(matches) <- c("Forecast", "Observed")
    out$matches <- matches
  
    X <- x$X.feats
    Y <- x$Y.feats

    X.feats <- Y.feats <- list()

    X.labeled <- Y.labeled <- matrix(0, xdim[1], xdim[2])

    if(verbose) cat("Loop through ", nmatches, " components of merges list to set up new (matched) features.\n")

    for(i in 1:nmatches) {

	if(verbose) cat(i, " ")
	tmp <- m[[ i ]]

	uX <- sort(unique(tmp[, 2]))
	uY <- sort(unique(tmp[, 1]))

	nX <- length(uX)
	nY <- length(uY)

	Xtmp <- X[[ uX[1] ]]
        Ytmp <- Y[[ uY[1] ]]

	if(nX > 1) for(j in 2:nX) Xtmp <- union.owin(Xtmp, X[[ uX[j] ]])
	if(nY > 1) for(k in 2:nY) Ytmp <- union.owin(Ytmp, Y[[ uY[k] ]])

	X.feats[[ i ]] <- Xtmp
	Y.feats[[ i ]] <- Ytmp

	X.labeled[ Xtmp$m ] <- i
	Y.labeled[ Ytmp$m ] <- i

	Xtmp <- Ytmp <- NULL

    } # end of for 'i' loop.

    unX <- sort(x$unmatched$X)
    unY <- sort(x$unmatched$Xhat)

    nX2 <- length(unX)
    nY2 <- length(unY)

    if(nX2 > 0) {

	if(verbose) cat("\nLoop to add/re-label all unmatched observed features.\n")

	vxunmatched <- (nmatches + 1):(nmatches + nX2)

	for(i in 1:nX2) {

	    Xtmp <- X[[ unX[i] ]]
	    X.feats[[ nmatches + i ]] <- Xtmp
	    X.labeled[ Xtmp$m ] <- nmatches + i

	} # end of for 'i' loop.

    } else vxunmatched <- integer(0)

    if(nY2 > 0) {

	if(verbose) cat("\nLoop to add/re-label all unmatched forecast features.\n")

        fcunmatched <- (nmatches + 1):(nmatches + nY2)

        for(i in 1:nY2) {

            Ytmp <- Y[[ unY[i] ]]
            Y.feats[[ nmatches + i ]] <- Ytmp
            Y.labeled[ Ytmp$m ] <- nmatches + i

        } # end of for 'i' loop.

    } else fcunmatched <- integer(0)

    out$X.feats <- X.feats
    out$Y.feats <- Y.feats
    out$X.labeled <- X.labeled
    out$Y.labeled <- Y.labeled
    out$unmatched <- list(X = vxunmatched, Xhat = fcunmatched)

    # 'out' should already have class "matched".
    # class(out) <- "matched"
    return(out)

} # end of 'MergeForce' function.
