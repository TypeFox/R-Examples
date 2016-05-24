KernSur <- function(
		 x, 
		 y, 
		 xgridsize=100, 
		 ygridsize=100,	 
		 correlation, 
		 xbandwidth, 
		 ybandwidth, 
		 range.x,
		 range.y,
		 na.rm=FALSE
		 ) 

{
# multiplier for the xbandwidth which gives the tails at either end for the 
# final correlated kernel density estimate only operates as default follows
# Wand and Jones' 1.5 h default
deadspace <- 1.5        
# number of cases to deal with
cases <- length(x)		
if(cases != length(y)){stop("x and y vectors of unequal lengths")}

# set up some control variables
flag1 <- 1; flag2 <- 1; flag3 <- 1; flag4 <- 1; flag5 <- 0

# vector of correlation coefficients
	 # default behaviour
         if(missing(correlation))
                {
		flag5 <- 0
		correlation <- correlationselect(x, y, correlation=FALSE, cases)
		} 
	if(flag5){correlation <- correlationselect(x, y, correlation, cases)}


# bandwidth selection
	# default bandwidths
        if(missing(xbandwidth))
                {
		flag1 <- 0
		xbandwidth <- bandwidthselect(x, bandwidths=FALSE)
		} 
        if(missing(ybandwidth))
                {
		flag2 <- 0
		ybandwidth <- bandwidthselect(y, bandwidths=FALSE)
		} 
	# user supplied bandwidths
	if(flag1){xbandwidth <- bandwidthselect(x, xbandwidth)}
	if(flag2){ybandwidth <- bandwidthselect(y, ybandwidth)}


# x-y related stuff done now do the NA handling
# put it all together into a data frame or na.omit doesn't work
z <- data.frame(x,y, correlation, xbandwidth, ybandwidth)
	# if NAs not allowed fail the function
	if(na.rm == FALSE){na.fail(z)}
	# get rid of NA cases
	if(na.rm == TRUE){z <- na.omit(z)}
# reassign the vectors with NAs removed
x <- z$x; y <- z$y; correlation <- z$correlation; xbandwidth <- z$xbandwidth; ybandwidth <- z$ybandwidth


# range selection and ordinate generation
	# default range of values
        if(missing(range.x))
                {
		flag3 <- 0
		xvals <- rangeselect(x, rnge=FALSE, xgridsize, xbandwidth, deadspace)
		}
       if(missing(range.y))
                {
		flag4 <- 0
		yvals <- rangeselect(y, rnge=FALSE, ygridsize, ybandwidth, deadspace)
		}
	# user supplied ranges
	if(flag3){xvals <- rangeselect(x, range.x, xgridsize, xbandwidth, deadspace)}
	if(flag4){yvals <- rangeselect(y, range.y, ygridsize, ybandwidth, deadspace)}


# setup ordinate vector lengths
ordsinx <- length(xvals)
ordsiny <- length(yvals)
# reset cases with NAs removed
cases <- length(x)
# generate the vector of squared correlations
correlationsq <- correlation ^ 2 
# define the  kernelsurface array
corker <- rep(0, (ordsinx * ordsiny)) 

# invoke the .c module 
out <- .C( 
	"GenKernSur",
	 as.double(corker), 
	 as.integer(ordsinx), 
	 as.integer(ordsiny),
	 as.double(x), 
	 as.double(y), 
	 as.double(xvals), 
	 as.double(yvals), 
	 as.double(xbandwidth), 
	 as.double(ybandwidth), 
	 as.double(correlation), 
	 as.double(correlationsq),
	 as.integer(cases),
	 PACKAGE="GenKern"
	 ) 
 
# assign the return values 
zden <- out[[1]] 
dim(zden) <- c(ordsinx, ordsiny) 

# condtruct list to return or R>1.8 complains
op <- list(xvals, yvals, zden)
names(op) <- c("xords", "yords", "zden")

return(op) 
}
