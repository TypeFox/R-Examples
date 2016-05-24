KernSec <- function(
		    x,
		    xgridsize=100,
		    xbandwidth,
		    range.x
		   )
{

deadspace <- 1.5
flag.type.II <- 0
flag.range <- 0


# let there be n observations, and xgridsize ordinates
# xbandwidth can be missing, 1, or n, or xgridsize
# range.x can be missing, 2, or xgridsize

# uniform bandwidth automatically selected
if(missing(xbandwidth))
	{
	xbandwidths <- bandwidthselect(as.vector(na.omit(x)), bandwidths=FALSE)

	if(missing(range.x))
		{
		xords <- rangeselect(x, rnge=FALSE, as.integer(xgridsize), xbandwidths, deadspace)
		if(any(is.na(x)) == TRUE){warning("NAs in x - removing: "); x <- as.vector(na.omit(x))}
		flag.range <- 1
		}


	if((! flag.range) && (! missing(range.x)) && (length(range.x) == 2))
		{
		xords <- rangeselect(x, rnge=range.x, as.integer(xgridsize), xbandwidths, deadspace)
		if(any(is.na(x)) == TRUE){warning("NAs in x - removing: "); x <- as.vector(na.omit(x))}
		flag.range <- 1
		}


	if((! flag.range) && (! missing(range.x)) && (length(range.x) > 5))
		{
		xords <- range.x
		if(any(is.na(xords)) == TRUE){warning("NAs in range.x - removing: "); x <- as.vector(na.omit(xords))}
		if(any(is.na(x)) == TRUE){warning("NAs in x - removing: "); x <- as.vector(na.omit(x))}
		flag.range <- 1
		}
	}



# uniform bandwidth passed by user
if((! missing(xbandwidth)) && (length(xbandwidth) == 1))
	{
	xbandwidths <- bandwidthselect(as.vector(na.omit(x)), bandwidths=xbandwidth)

	if(missing(range.x))
		{
		xords <- rangeselect(x, rnge=FALSE, as.integer(xgridsize), xbandwidths, deadspace)
		if(any(is.na(x)) == TRUE){warning("NAs in x - removing: "); x <- as.vector(na.omit(x))}
		flag.range <- 1
		}


	if((! flag.range) && (! missing(range.x)) && (length(range.x) == 2))
		{
		xords <- rangeselect(x, rnge=range.x, as.integer(xgridsize), xbandwidths, deadspace)
		if(any(is.na(x)) == TRUE){warning("NAs in x - removing: "); x <- as.vector(na.omit(x))}
		flag.range <- 1
		}


	if((! flag.range) && (! missing(range.x)) && (length(range.x) > 5))
		{
		xords <- range.x
		if(any(is.na(xords)) == TRUE){warning("NAs in range.x - removing: "); x <- as.vector(na.omit(xords))}
		if(any(is.na(x)) == TRUE){warning("NAs in x - removing: "); x <- as.vector(na.omit(x))}
		flag.range <- 1
		}
	}





# variable bandwidth passed by user one for each observation
if((! missing(xbandwidth)) && (length(xbandwidth) == length(x)))
	{
	xbandwidths <- xbandwidth

	z <- data.frame(x, xbandwidth)
	if(any(is.na(z)) == TRUE){warning("NAs in bandwidths or/and data - removing: "); z <- na.omit(z)}
	x <- z$x; xbandwidths <- z$xbandwidth

	if(missing(range.x))
		{
		xords <- rangeselect(x, rnge=FALSE, as.integer(xgridsize), xbandwidths, deadspace)
		flag.range <- 1
		}


	if((! flag.range) && (! missing(range.x)) && (length(range.x) == 2))
		{
		xords <- rangeselect(x, rnge=range.x, as.integer(xgridsize), xbandwidths, deadspace)
		flag.range <- 1
		}


	if((! flag.range) && (! missing(range.x)) && (length(range.x) > 5))
		{
		xords <- range.x
		if(any(is.na(xords)) == TRUE){warning("NAs in range.x - removing: "); x <- as.vector(na.omit(xords))}
		flag.range <- 1
		}
	}



# variable bandwidth one for each ordinate
if((! missing(xbandwidth)) && (! missing(range.x)) && (length(xbandwidth) == length(range.x)))
	{
	flag.type.II <- 1
	z <- data.frame(range.x, xbandwidth)
	if(any(is.na(z)) == TRUE){warning("NAs in bandwidths or/and ordinates - removing: "); z <- na.omit(z)}
	xords <- z$range.x; xbandwidths <- z$xbandwidth
	if(any(is.na(x)) == TRUE){warning("NAs in x - removing: "); x <- as.vector(na.omit(x))}
	flag.range <- 1
	}



# variable bandwidth one for each ordinate
if((! missing(xbandwidth)) && (length(xbandwidth) == xgridsize))
	{
	flag.type.II <- 1
	if(any(is.na(x)) == TRUE){warning("NAs in x - removing: "); x <- as.vector(na.omit(x))}
	xbandwidths <- xbandwidth
	if(any(is.na(xbandwidths)) == TRUE){warning("NAs in bandwidths - removing: "); xbandwidths <- as.vector(na.omit(xbandwidth))}


	if(missing(range.x))
		{
		xords <- rangeselect(x, rnge=FALSE, as.integer(xgridsize), xbandwidths, deadspace)
		flag.range <- 1
		}


	if((! flag.range) && (! missing(range.x)) && (length(range.x) == 2))
		{
		xords <- rangeselect(x, rnge=range.x, as.integer(xgridsize), xbandwidths, deadspace)
		flag.range <- 1
		}

	}


# trap anything else
if((! missing(xbandwidth)) && (! missing(range.x)))
#{if(length(x) == length(range.x) == length(xbandwidth)){stop("too many inputs of same length")}}
{if((length(x) == length(range.x)) && (length(x) == length(xbandwidth))){stop("too many inputs of same length")}}

if(! flag.range){stop("undefined problem with input")}

# generate a vector of length xordslen with zeros in it to contain
# the density estimate
xordslen <- length(xords)
est <- rep(0, xordslen)

# invoke the .c module
out <- .C(
	 "GenKernSec",
	 as.double(x),
	 as.integer(length(x)),
	 as.double(xords),
	 as.double(xbandwidths),
	 as.double(est),
	 as.integer(xordslen),
	 as.integer(flag.type.II),
	 PACKAGE="GenKern"
	 )

# assign the return values
yden <- out[[5]]

# return list as from R 1.8 or awkward messages from the R 
# interpreter - don't forget the list names
op <- list(xords, yden); names(op) <- c("xords", "yden")

return(op)
}

