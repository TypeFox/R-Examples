# function nearest which picks the index of a vector which is nearest to a value
# the vector doesn't have to be in any particular order this routine will just
# give the nearest numbers index
# you can specify whether to accept values which are outside the range of the vector
# the default is not to but if you do then you can be way outside the range and
# get a return value
# the only bug is that if you don't specify that values which are outside the
# range of the vector are acceptable then any value which is equal to the minimum
# or maximum of the vector will also return an error so if too many errors
# are returned call with the outside=T flag enabled
# this function is fairly bullet proof but don't expect it to get things right when
# there are many tied values in the vector from playing with it the function returns
# vector of the tied values - so be prepared

nearest <- function(x, xval, outside=FALSE, na.rm=FALSE)
{

# Do the NA handling
# put it all together into a data frame or na.omit doesn't work
x1 <- x
z <- data.frame(x, x)
# if NAs not allowed fail the function
if(na.rm == FALSE){na.fail(z)}
# get rid of NA cases
if(na.rm == TRUE){z <- na.omit(z)}
# reassign the vectors with NAs removed
x <- z$x


# if the value is outside the range of the vector and it isn't acceptable
# then issue an error and stop
if(outside == FALSE){if((max(x) <= xval) || (min(x) >= xval)) {stop("value outside vector range")}}

# if the value is outside the range of the vector and this is acceptable then
# merely assign one of the index with one of the extreme values in it

if(outside == TRUE)
	{
	if((max(x) <= xval) || (min(x) >= xval))
		{
		sorx <- sort(x)
		if(abs(sorx[1] - xval) < abs(sorx[length(sorx)]- xval))
			{index <- 1; vally <- sorx[index]; index <- which(x1 == vally); return(index)}

		if(abs(sorx[1] - xval) > abs(sorx[length(sorx)]- xval))
			{index <- length(sorx); vally <- sorx[index]; index <- which(x1 == vally); return(index)}
		}
	}


# for most cases in which the value falls within the vector find the nearest
# value and assign that index to the return value
sorx <- sort(x)
upp <- which(sorx >= xval)[1]
low <- which(sorx <= xval); low <- low[length(low)]
upp <- sorx[upp]; low <- sorx[low]

if(upp == low) {index <- which(x == upp)}
if((abs(upp - xval)) >= (abs(low - xval))) {index <- which(x1 == low)}
if((abs(upp - xval)) < (abs(low - xval))) {index <- which(x1 == upp)}

return(index)
}


