hwt <-
function (x, type=c("wavelet", "station"), reindex=FALSE)
{

    type <- match.arg(type)

    n <- length(x)

    if (n==0)
	stop("Vector x has to have some entries.")

    else if (any(is.nan(x)))
	stop("All entries in x have to be numbers")

    ans <- NULL

    nsteps <- floor(log(n, base=2))

    ansd <- ansc <- vector("list", nsteps)

    vc <- x

    h <- c(1, 1)/sqrt(2)
    g <- c(-1, 1)/sqrt(2)
 
    if (type=="wavelet")	{
	    for(i in 1:nsteps)	{

		vd <- as.vector(filter(x=vc, filter=g))
		vc <- as.vector(filter(x=vc, filter=h))

		vc <- vc[seq(from=1, to=length(vc), by=2)]
		vd <- vd[seq(from=1, to=length(vd), by=2)]

		ansc[[nsteps-i+1]] <- vc[!is.na(vc)]
		ansd[[nsteps-i+1]] <- vd[!is.na(vd)]
		}
	}
    else if (type=="station")	{

	    work.h <- h
	    work.g <- g

	    for(i in 1:nsteps)	{

		ansc[[nsteps-i+1]] <- as.vector(filter(x=vc, filter=work.h))
		work.h <- as.vector(filter(x=zeropad(work.h), filter=h, circular=TRUE, sides=1))

		ansd[[nsteps-i+1]] <- as.vector(filter(x=vc, filter=work.g))
		work.g <- as.vector(filter(x=zeropad(work.g), filter=h, circular=TRUE, sides=1))

		}
	}

#
# Since this object is the result of HWT on a non-dyadic number of observations
# one can compare it to an object which is the dyadic HWT of the vector which
# is the next high dyadic vector in length. E.g. if the length here is 30
# then the next highest will be 32. Sometimes we want this object to behave
# like the 32 object, and so we have to insert an extra coarser scale, but
# blank. We do this if reindex=TRUE 
#
  if (reindex==TRUE)	{

    ransd <- ransc <- vector("list", nsteps+1)

    ransd[1] <- list(NULL)
    ransc[1] <- list(NULL)

    for(i in 1:nsteps)	{

	ransd[[i+1]] <- ansd[[i]]
	ransc[[i+1]] <- ansc[[i]]

	}
    ansd <- ransd
    ansc <- ransc
    nsteps <- nsteps+1
    }


#
# For ansc and ansd level 1 corresponds to the coarsest scale and level
# nsteps to the finest. This is like wd in wavethresh but the scale here
# is one more than in wavethresh (as wavethresh indexes scales from 0=coarsest,
# not 1=coarsest.
#
l <- list(c=ansc, d=ansd, nlevels=nsteps, type=type, reindex=reindex)

class(l) <- "hwtANYN"

return(l)

	
}
