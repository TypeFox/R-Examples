#--------------------------------------------------------------------
#   locpol.R (npsp package demo)
#--------------------------------------------------------------------
# EXAMPLES:
#   locpol.bin  S3 class and methods
#   locpol()    S3 generic
#
#   (c) R. Fernandez-Casal         Last revision: Mar 2014
#--------------------------------------------------------------------
library(npsp)

#--------------------------------------------------------------------
# Regularly spaced 2D data
#--------------------------------------------------------------------
nx <- c(40, 40)   # ndata = prod(nx)
x1 <- seq(-1, 1, length.out = nx[1])
x2 <- seq(-1, 1, length.out = nx[2])
trend <- outer(x1, x2, function(x,y) x^2 - y^2)
spersp( x1, x2, trend, main = 'Trend', zlab = 'y')

set.seed(1)
y <- trend + rnorm(prod(nx), 0, 0.1)
x <- as.matrix(expand.grid(x1 = x1, x2 = x2)) # two-dimensional grid
spersp( x1, x2, y, main = 'Data')

#  # Binning
#  bin <- binning(x, y)
#  str(bin)
#  dim(bin)
#  dimnames(bin)
#  str(coords(bin))
#  str(coordvalues(bin))
#  simage(bin, main = 'Binning surface')


# Local polynomial kernel regression 
lp <- locpol(x, y, h = diag(c(0.3, 0.3)))
str(lp)
dim(lp)
dimnames(lp)
str(coords(lp))
str(coordvalues(lp))

spersp(lp, main = 'locpol surface', zlab = 'y')

# lopol from binned data
lp2 <- locpol(lp, h = diag(c(0.1, 0.1)))  # avoids binning
spersp(lp2, main = 'locpol surface (h=0.1)', zlab = 'y')


#--------------------------------------------------------------------
# Regularly spaced 1D data
#--------------------------------------------------------------------
# one-dimensional data grid
nx <- 1000
x <- seq(0, 1, length.out = nx)
f1d <- function(x) sin( 2*pi*x )
y <- f1d(x) + rnorm(nx, 0, 0.5)

plot(x, y, type='l', col='darkgray')
curve(f1d, 0, 1, add = TRUE)

lp <- locpol(x, y, h = 0.05, nbin = 100)
lines(coords(lp), lp$biny, lwd = 2, lty = 2)
lines(coords(lp), lp$est, lwd = 2)

str(lp)

# cuadratic fit + derivatives

lp <- locpol(x, y, h = 0.3, nbin = 100, drv = TRUE)
lines(coords(lp), lp$est, lwd = 2, col = 2)

plot(coords(lp), lp$deriv, type = "l", main = "Estimated first derivative")