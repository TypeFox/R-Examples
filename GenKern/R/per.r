# per.r find the value of the ith percentage point of a probability distribution
# by assuming that if the vals has the same number of elements as the density then
# these must correspond to an approximation for bin centres and interpolates by
# splitting the difference between the centres to find bin extremes
# if the vals vector is and element greater than the density vector then it assumes
# the vals vector represents an approximation to the bin extremes so leaves them alone
# NAs are allowed in the density - if na.rm is TRUE then it assumes that they're really
# zero values - if na.rm is FALSE then it fails with warning
# NAs are not allowed in the vals vector at all - plus there can be no repeated values
# in it, and they MUST be in ascending order - contravene any of these and the function
# will fail with warning
# normally the density cannot have negative values, but if neg.rm is set to TRUE
# then all negative values from the density will be assumed to be zero - this
# occasionally happens from various kde functions

per <- function(den, vals, point, na.rm=FALSE, neg.rm=FALSE)
{
# handle various NA actions
# fail if they appear at all
if(na.rm == FALSE){na.fail(den); na.fail(vals)}
# assume NA values in a pdf are really zero - fail NAs from ordinates
if(na.rm == TRUE){den[which(den == "NA")] <- 0; na.fail(vals)}

lenden <- length(den)
lenvals <- length(vals)

# check we have the vectors of the right lengths approximate to bin centres or bin limits
if((lenvals < lenden) || (lenvals > (lenden + 1))){stop("mismatch in vector lengths")}
# test to see whether there are any repeated ordinate values
if(length(unique(vals)) != lenvals){stop("repeated values in ordinates")}
# make sure the ordinates are arranged in ascending order
if(length(unique(vals == sort(vals))) == 2){stop("ordinate values not in ascending order")}

# negative density values
# fail if there are any density values less than zero
if(neg.rm == FALSE){if(any(den < 0)){stop("negative values in the density vector")}}
# convert all negative values of the density to zero if required behaviour
if(neg.rm == TRUE){den[which(den < 0)] <- 0}


# if vals is a vector approximating bin centres calculate the bin extremes
	if(lenden == lenvals)
		{
		extremes <- rep(0, lenden + 1)
		outty <- .C("getlims", as.double(vals), as.double(extremes), as.integer(lenden), PACKAGE="GenKern")
		# reassign the bin extreme vector
		extremes <- outty[[2]]
		}

# if vals is a vector approximating bin extremes then set up the extremes vector
if(lenvals == (lenden + 1)){extremes <- vals}
		
#####################################################################################

lenextremes <- length(extremes)

# calculate vectors of extreme values so a diffrence relating to each
# value of the density may be calculated
extremes1 <- extremes[1:(lenextremes - 1)]
extremes2 <- extremes[2:lenextremes]
# use the difference to calculate the integrated area for each bin
volumes <- (abs(extremes2 - extremes1)) * den
# can now legitimately normalise the volumes
sumvolumes <- sum(volumes)
volumes <- volumes / sumvolumes

# find the first bin for which the integrated area exceeds that of the
# requested point
culvolumes <- cumsum(volumes)
bin <- which(culvolumes > point)[1]

# calculate the distance across the bin which the integrated area
# equals the point which specified by point
# notice special case for if the first bin is the one in which
# the point value occurs
if(bin == 1){prop <- point / culvolumes[bin]}
if(bin > 1){prop <- (point - culvolumes[bin-1]) / (culvolumes[bin] - culvolumes[bin -1])}

# prop cannot be greater then 1 - if it is then the algorthm is duff
if(abs(prop) > 1){stop("Bugger! - duff code again - notify developers")}

# work out the interpolated value as prop * distance between the two bin
# extremes
xvalue <- extremes[bin] + ((extremes[bin + 1] - extremes[bin]) * prop)
return(xvalue)
}

