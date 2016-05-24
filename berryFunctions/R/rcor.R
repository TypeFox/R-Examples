#' Spatially correlated random values
#' 
#' Generate random values, but with spatial correlation
#' 
#' @author Berry Boesenkool, \email{berry-b@@gmx.de}, Jan 2016
#' @references For regular grid fields, see: \url{http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/}
#' @note This function is not yet running
#' @export
#' @examples
#' # ToDO: expand function + write examples
#'
#' @param n Number of randomly scattered points to be created
#' @param xmin  Smallest x-coordinate
#' @param xmax Largest x-coordinate
#' @param ymin Smallest y-coordinate
#' @param ymax Largest y-coordinate
#' @param allplots Save all plots in dummy.pdf? DEFAULT: TRUE
#' @param zstart Vector of starting values for seeding some points, must be shorter than n
#' @param zfun A function for random noise creation, taking n as first argument. DEFAULT: rnorm
#' @param \dots Further arguments passed to zfun, like sd=3
#'
rcor <- function(
n,
xmin,
xmax,
ymin,
ymax,
allplots=TRUE,
zstart,
zfun=rnorm,
...)
{
if(allplots) { pdf("dummy.pdf") ; on.exit(dev.off())}
# coordinate vectors
x <- runif(n, xmin, xmax)
y <- runif(n, ymin, ymax)
z <- rep(NA, n)
# number of seeded points:
ns <- length(zstart)
z[1:ns] <- zstart
# browser()
#### distances for each seeded point:
# compute autocorrelated z values:
for(i in (ns+1):n)
{
if(allplots) colPoints(x,y,z, add=FALSE, asp=1)
# get value from closest populated point and add noise
d <- distance(x[i],y[i], x,y)
d[is.na(z)] <- NA
z[i] <- z[which.min(d)] + zfun(n=1, ...)
}
if(allplots) colPoints(x,y,z, add=FALSE, asp=1)
# function output
data.frame(x,y,z)
}

if(FALSE){
coor <- rcor(n=100, xmin=45, xmax=80, ymin=230, ymax=270, zstart=c(20,30,40), sd=3)
colPoints(x,y,z, data=coor, add=FALSE, asp=1, col=divPal(100))
points(coor$x[1:3], coor$y[1:3], cex=1.4)
}

