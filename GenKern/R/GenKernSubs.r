# subroutines for KernSec and KernSur
# written because much of the checking between the two functions
# was common and it seemed so much more elegant to have all data
# verification functions in a single file



# generates a suitable vector of bandwidths for dat
# from any allowed input from bandwidths which can be
# 1 - a vector of length equal to that of dat
# 2 - a single value
# 3 - nothing - if true bandwidths should be set to FALSE
bandwidthselect <- function(dat, bandwidths)
{

cases <- length(dat)
# calculate the default bandwidths
if(any(is.na(bandwidths))){bandwidths == FALSE}

if(! is.numeric(bandwidths))
	{
	# dummy vector to get na.omit to work set up dat2 as repeated unit
	dat1 <- dat
	dat2 <- na.omit(data.frame(dat, dat1))$dat
	# trap vectors of x likely to upset dpik()
	if(length(unique(dat2)) < 5){stop("too few unique values - no auto bandwidth possible")}
	# calculate a vector of bandwidths based on dat using dpik()
	bandwidths <- rep(dpik(dat2), cases)
	}

# user supplied vector of bandwidths do nothing
bands <- length(bandwidths)

# user supplied single value bandwidths
if(bands == 1){bandwidths <- rep(bandwidths, cases)}

# trap anything unexpected
if(length(bandwidths) != cases){stop("bandwidth vector of funny size")}
if(any(bandwidths == "NA")){stop("bandwidth vector contains NA's")}
return(bandwidths)
}









# outputs a vector of ordinates at which to calculate the density
# users may specify the range of the extreme ordinates, give
# a vector of ordinates themselves, or nothing, in which the
# default is gridsize equally spaced points with deadspace times the
# maximum bandwidth at each end from the data
rangeselect <- function(dat, rnge, gridsize, bandwidth, deadspace)
{

# default behaviour
if(rnge[1] == "FALSE")
	{
	max.dat <- max(dat, na.rm=TRUE)
	min.dat <- min(dat, na.rm=TRUE)
	# make sure that all the values are not the same
	if((max.dat - min.dat) == 0){stop("all x's are the same - no auto range possible")}
	maxbw <- max(bandwidth, na.rm=TRUE)
	# calculate a range based upon deadspace and bandwidths
	rnge <- c( (min.dat - (deadspace * maxbw)),(max.dat + (deadspace * maxbw)) )
	# calculate the vector of ordinates
	ords <- seq(rnge[1], rnge[2], length=round(gridsize))
	}



if(rnge[1] != "FALSE")
{
# if user has sent values for the extreme limits of the range
if(length(rnge) == 2)
	{
	rnge <- sort(rnge)
	ords <- seq(rnge[1], rnge[2], length=as.integer(gridsize))
	}

# if user has specified their own ordinates
if(length(rnge) > 2){ords <- sort(rnge)}
}


# test to see whether there are any repeated ordinate values from a user vector
if(length(unique(ords)) != length(ords)){stop("repeated values in ordinates")}

return(ords)
}



# constructs a vector of correlations to use with KernSur
# users can:
# 1 - send a single correlation of their choise
# 2 - a vector of correlations for each case
# 3 - or nothing in which case a default correlation is established
correlationselect <- function(x, y, correlation, cases)
{


# default if correlation isn't there
if(correlation == "FALSE")
	{
	min.x <- min(x, na.rm=TRUE)
	max.x <- max(x, na.rm=TRUE)
	min.y <- min(y, na.rm=TRUE)
	max.y <- max(y, na.rm=TRUE)

	if((min.x - max.x == 0) || (min.y - max.y == 0))
		{corre <- 0}
	
	if((min.x - max.x != 0) && (min.y - max.y != 0))
		{corre <- cor(x,y, use="complete.obs")}
		
	correlation <- rep(corre, cases)
	}

# if user sends a single value
if(length(correlation) == 1)
	{
	corre <- correlation
	correlation <- rep(corre, cases)
	} 

# reset untenable values for the correlation
if(length(correlation) == cases)
	{
	correlation[which(correlation > 0.999)] <- 0.999
	correlation[which(correlation < -0.999)] <- -0.999
	}

# trap bad vectors for the correlation
if(length(correlation) != cases){stop("wrong vector length for correlation")}

return(correlation)
}



# getlims.r return a vector of bin limits based on a vector of bin centres
getlims <- function(centres)
{
# set up the vector length
bins <- length(centres)

# test to see whether there are any repeated bin centre values
if(length(unique(centres)) != bins){stop("repeated values in bin centres")}

# make sure the bin centres are arranged in ascending order
if(length(unique(centres == sort(centres))) == 2){stop("bin centre values not in ascending order")}

# define a vector of bin limits to be a vector with n+1 bins
limits <- rep(0, bins+1)

# calculate the bin limits
outty <- .C("getlims", as.double(centres), as.double(limits), as.integer(bins), PACKAGE="GenKern")

# reassign the output
limits <- outty[[2]]
return(limits)
}








