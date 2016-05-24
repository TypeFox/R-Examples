#$Id: lshorth.R 49 2008-02-11 12:39:48Z gsawitzki $
lshorth <-
function (x, probs = NULL,  plot = TRUE, na.rm=FALSE, ...) 
{
    if (!is.numeric(x)) 
        stop("'x' must be numeric")
    xname <- paste(deparse(substitute(x), 50), collapse = "\n")
    # name must be set before removin NAs 

	if (na.rm) {x <- x[is.finite(x)]} else {
		if ( any(!is.finite(x))) stop("'x' contains infinite or missing values. Not yet supported.")}

    if (length(x)<1)
	    stop("'x' must have positive length")
	count <- length(x)
	# if (is.null(psteps)) {psteps <- 2*ceiling(log2(length(x))+1)}
	
		if (is.null(probs)){
			ppx <-ceiling(log2(count) / 2) 
			probs <- c( 1/2^(ppx:1),1-1/2^(2:ppx))
			} else {if (!is.numeric(probs)) stop("'probs' must be numeric")}
	# if (is.null(probs)) {probs <- (1:psteps)/(psteps + 1)}
    
    shorthm <- matrix(ncol = length(probs), nrow = length(x))
    xsort <- sort(x)
   for (px in 1:length(probs)) {
        {
            
            if (probs[px] <= 1/count) {
                shorthvec <- rep.int(0, count)
            }
            else {
                jDelta <- ceiling(probs[px] * count) - 1
                 jmax <- count - jDelta
                lenvec <- xsort[-(1:(jDelta))] - xsort[-((count - 
                  jDelta + 1):count)]
                shorthvec <- vector(mode = "numeric", count)
                minlenj <- 1
                shorthvec[1] <- minlen <- lenvec[1]
                for (i in 2:count) {
                  jmin <- i - jDelta
                  if (jmin < 1) {
                    jmin <- 1
                  }
                  if (i + jDelta > count) {
                    minlenj <- jmin - 1 + which.min(lenvec[jmin:(count - 
                      jDelta)])
                  }
                  else {
                    if (minlenj < jmin) {
                      minlenj <- jmin - 1 + which.min(lenvec[jmin:i])
                    }
                    else {
                      if (lenvec[i] < minlen) {
                        minlenj <- i
                      }
                    }
                  }
                  shorthvec[i] <- minlen <- lenvec[minlenj]
                }
            }
            shorthm[, px] <- shorthvec
        }
    }
    shorthm <- structure(list(x = xsort, lshorth = shorthm, 
        probs = probs, xname=xname), class = "lshorth")
    if (plot) {
        plot(shorthm, probs=probs,  ...)
    }
    invisible(shorthm)
}
