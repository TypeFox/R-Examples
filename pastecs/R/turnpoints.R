"turnpoints" <-
function(x) {
	data <- deparse(substitute(x))
	if (is.null(ncol(x)) == FALSE)
		stop("Only one series can be treated at a time")
    x <- as.vector(x)
    n <- length(x)
    diffs <- c(x[1]-1, x[1:(n-1)]) != x
    uniques <- x[diffs]
    n2 <- length(uniques)
    poss <- (1:n)[diffs]
    exaequos <- c(poss[2:n2], n+1) - poss - 1
    if (n2 < 3) {			#We need at least 3 unique values!!!
    	warning("Less than 3 unique values, no calculation!")
    	nturns <- NA
    	firstispeak <- FALSE
    	peaks <- rep(FALSE, n2)
    	pits <- rep(FALSE, n2)
    	tppos <- NA
    	proba <- NA
    	info <- NA
    } else {
    	# The following code is faster in R, but do not work all the time!
    	#	ex <- embed(uniques, 3)	# Works only in R!
    	#	peaks <- c(FALSE, max.col(ex) == 2, FALSE)
    	#	pits <- c(FALSE, max.col(-ex) == 2, FALSE)
    	m <- n2 - 2
    	ex <- matrix(uniques[1:m + rep(3:1, rep(m, 3)) - 1], m)
    	peaks <- c(FALSE, apply(ex, 1, max, na.rm=TRUE) == ex[, 2], FALSE)
    	pits <- c(FALSE, apply(ex, 1, min, na.rm=TRUE) == ex[, 2], FALSE)
    	tpts <- peaks | pits
    	if (sum(tpts) == 0) {	# No turning point
    		nturns <- 0
    		firstispeak <- FALSE
			peaks <- rep(FALSE, n2)
			pits <- rep(FALSE, n2)
			tppos <- NA
			proba <- NA
    	    info <- NA
    	} else {
    		tppos <- (poss + exaequos)[tpts] 	# This way, we consider the last element of duplicates, as in PASSTEC 2000
    		tptspos <- (1:n2)[tpts] 
    		firstispeak <- tptspos[1] == (1:n2)[peaks][1] 
    		nturns <- length(tptspos)
    		if (nturns < 2) {
    			inter <- n2 + 1
    			posinter1 <- tptspos[1]
    		} else {
    			inter <- c(tptspos[2:nturns], n2) - c(1, tptspos[1:(nturns-1)]) + 1
    			posinter1 <- tptspos - c(1, tptspos[1:(nturns-1)])
    		}
    		posinter2 <- inter - posinter1
			posinter <- pmax(posinter1, posinter2)
			proba <- 2 / (inter * gamma(posinter) * gamma(inter - posinter + 1))
    		info <- -log(proba, base=2)
    	}
    }
    res <- list(data=data, n=n, points=uniques, pos=(poss + exaequos), exaequos=exaequos, nturns=nturns, firstispeak= firstispeak, peaks=peaks, pits=pits, tppos=tppos, proba=proba, info=info)
    class(res) <- "turnpoints"
    res
}
