"isotone" <-
function (x, wt = rep(1, length(x)), increasing = FALSE, Ccode=TRUE) 
{
    nn <- length(x)

    if (Ccode==TRUE)	{
		ans <- .C("isotoneC",
			x = as.double(x),
			wt = as.double(wt),
			nn = as.integer(nn),
			increasing = as.logical(increasing),
			PACKAGE="DDHFm")
		return(ans$x)
    }
    else	{
	    if (nn == 1) 
		return(x)
	    if (!increasing) 
		x <- -x
	    ip <- (1:nn)
	    dx <- diff(x)
	    nx <- length(x)
	    while ((nx > 1) && (min(dx) < 0)) {

		#cat("diff: ", dx, "\n")
		jmax <- (1:nx)[c(dx <= 0, FALSE) & c(TRUE, dx > 0)]
		jmin <- (1:nx)[c(dx > 0, TRUE) & c(FALSE, dx <= 0)]
		#cat("jmax: ", jmax, "\n")
		#cat("jmin: ", jmin, "\n")
		for (jb in (1:length(jmax))) {
		    ind <- (jmax[jb]:jmin[jb])
		    wtn <- sum(wt[ind])
		    x[jmax[jb]] <- sum(wt[ind] * x[ind])/wtn
		    wt[jmax[jb]] <- wtn
		    x[(jmax[jb] + 1):jmin[jb]] <- NA
		}
		ind <- !is.na(x)
		x <- x[ind]
		wt <- wt[ind]
		ip <- ip[ind]
		dx <- diff(x)
		nx <- length(x)
	    }
	    #print(x)
	    #print(ip)
	    jj <- rep(0, nn)
	    jj[ip] <- 1
	    z <- x[cumsum(jj)]
	    if (!increasing) 
		z <- -z
	    return(z)
	}
}

