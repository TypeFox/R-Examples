mnnd <- function(x,y) {
#
# Function mnnd to calculate the mean nearest neighbour distance
# between the points whose coordinates are stored in x and y.
#

n <- length(x)
if(n!=length(y)) stop('data lengths do not match')
dmb <- (max(x)-min(x))**2 + (max(y)-min(y))**2

.Fortran(
	"mnnd",
	x=as.double(x),
	y=as.double(y),
	n=as.integer(n),
	dmb=as.double(dmb),
	d=double(1),
	PACKAGE='deldir'
	)$d
}
