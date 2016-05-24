AntiAR <-
function (S) 
{
#
# Function takes spectral object S and halves its width (as it was computed
# on a data set that was reflected in its RH end to mitigate boundary effects.
#
# Get the number of levels, and put father and mother coefficients into matrices
#
J <- S$nlevels
Cm <- matrix(S$C, nrow = J + 1, byrow = TRUE)
Dm <- matrix(S$D, nrow = J, byrow = TRUE)
#
# Remove the right-hand side of the matrices
#
Cm <- Cm[1:J, 1:(ncol(Cm)/2)]
Dm <- Dm[1:(J - 1), 1:(ncol(Dm)/2)]
#
# Adjust the first-last database to take account of the new dimensions of
# the returned vectors.
#
fl.dbase <- first.last(LengthH = length(S$filter$H),
	DataLength = 2^(J - 1), type = "station", bc = "periodic",
	current.scale = 0)
#
# Make return object to be of class wd with the reduced coefficients
#
l <- list(C = as.vector(t(Cm)), D = as.vector(t(Dm)), nlevels = J - 
	1, fl.dbase = fl.dbase, filter = S$filter, type = "station", 
	bc = "periodic", date = S$date)
#
# Ensure returned object has correct class
#
class(l) <- "wd"
#
# Return object
#
return(l)
}
