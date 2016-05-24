# functions to assist in the creation of featuresets for blighty()
# Copyright - David Lucy January 2006

# converts fig files with the preamble removed to long
# lists of x-y coordinate pairs
#
# basically get your .fig file from xfig which looks like
###################################################################
# #FIG 3.2
# Landscape
# Center
# Metric
# A4
# 100.00
# Single
# -2
# 1200 2
# 2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 73
# 	 3440 5190 3459 5171 3459 5141 3459 5115 3459 5085 3482 5085
# 	 3508 5085 3523 5096 3545 5089 3549 5066 3568 5059 3587 5040
# 	 3598 5025 3613 5010 3628 4999 3643 4980 3662 4976 3680 4961
# 	 3703 4950 3725 4920
###################################################################
# and top and tail it so it looks like
###################################################################
# 	 3440 5190 3459 5171 3459 5141 3459 5115 3459 5085 3482 5085
# 	 3508 5085 3523 5096 3545 5089 3549 5066 3568 5059 3587 5040
# 	 3598 5025 3613 5010 3628 4999 3643 4980 3662 4976 3680 4961
# 	 3703 4950 3725 4920
###################################################################
# then complete the final row of coordinates with NA's so it is
###################################################################
# 	 3440 5190 3459 5171 3459 5141 3459 5115 3459 5085 3482 5085
# 	 3508 5085 3523 5096 3545 5089 3549 5066 3568 5059 3587 5040
# 	 3598 5025 3613 5010 3628 4999 3643 4980 3662 4976 3680 4961
# 	 3703 4950 3725 4920 NA   NA   NA   NA   NA   NA   NA   NA
###################################################################
# then this function can be applied to the now transformed .fig file
# you will need to supply the filename and assign it to a matrix
getdata <- function(filename)
{
# open the file and read in the data from a fig file with
# the top removed and any incomplete lines filled with NA's
dat <- read.table(filename, na="NA")

# set up a vector of the first file of the fig file
vec <- as.vector(dat[1,], mode="numeric")

# add to the vector all the other elements in the fig file
if(nrow(dat) >= 2)
	{
	for(ctr in 2:nrow(dat))
		{vec <- c(vec, as.vector(dat[ctr,], mode="numeric"))}

	}

# remove all NA's from the vector
vec <- vec[which(vec != "NA")]

# the vector of coordinates is made from the vector
coords <- matrix(vec, ncol=2, byrow=TRUE, dimnames=NULL)
colnames(coords) <- c("x","y")

return(coords)
}





invertfeature <- function(feature, allpoints)
{
# This function inverts the y's you get from xfig and from getdata()
# xfigs y coords are from top to bottom so we merely subtract each y coord
# for the feature set from the maximum y coord for all the feature sets
# under consideration

ymax <- max(allpoints[,2])
feature[,2] <- (ymax - feature[,2])

return(feature)
}




rescale <- function(feature, allpoints, limits, featuretype)
{
# This function rescales the vectors from invertfeature() to fit the
# extreme limits of the Cartesian coordinate system they're supposed to fill
# feature is a matrix of x-y coords from the invertfeature() function above
# allpoints is an rbinded matrix of all the features in the set
# limits are the maximum x and y limits of the featureset in allpoints
# in the final coordinate system
# can be limits <- c(min(x), max(x), min(y), max(y)) although I tend to
# take them from the coordinate system prescribed by the ordinance survey

x1 <- min(allpoints[,1])
x2 <- max(allpoints[,1])
y1 <- limits[1]
y2 <- limits[2]
feature[,1] <- (((feature[,1] - x1) / (x2 - x1)) * (y2 - y1)) + y1

x1 <- min(allpoints[,2])
x2 <- max(allpoints[,2])
y1 <- limits[3]
y2 <- limits[4]
feature[,2] <- (((feature[,2] - x1) / (x2 - x1)) * (y2 - y1)) + y1

# put in the feature type designation as the first x data point
# 1 is an area
# 2 is an inverted area
# 3 is a linear feature
feature <- rbind(c(featuretype, 0), feature)

return(feature)
}


