#--------------------------------------------------------------------
#   binning.R (npsp package demo)
#--------------------------------------------------------------------
# EXAMPLES:
#   binning   S3 generic
#   bin.data  S3 class and methods
#   locpol    S3 generic  (locpol.bin.data)
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

# Binning
bin <- binning(x, y)
str(bin)
dim(bin)
dimnames(bin)
str(coords(bin))
str(coordvalues(bin))
simage(bin, main = 'Binning surface')

# Local polynomial kernel regression
lp <- locpol(bin, h = diag(c(0.3, 0.3)), hat.bin = TRUE)
## Equivalent to:  lp <- locpol(x, y, h = diag(c(0.3, 0.3)) # not equal...
str(lp)
spersp(lp, main = 'locpol surface', zlab = 'y')

# Estimation with (binning) hat matrix
if (!is.null(lp$locpol$hat)) {
    y2 <- trend + rnorm(prod(nx), 0, 0.1)  # new y-data
    bin2 <- binning(x, y2)
    est2 <- lp$locpol$hat %*% as.numeric(bin$biny)
    est2 <- data.grid(est2, grid = bin2$grid)
    spersp(est2, main = 'locpol surface (binning hat matrix)', zlab = 'y')
} else cat("'locpol.bin' object must be created with parameter 'hat.bin = TRUE'")
