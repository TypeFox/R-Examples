#--------------------------------------------------------------------
#   datagrid.R (npsp package demo)
#--------------------------------------------------------------------
# EXAMPLES:
#   data.grid   S3 class and methods
#   grid.par    S3 class and methods
#   coords      S3 generic
#   coordvalues S3 generic
#
#   (c) R. Fernandez-Casal         Last revision: Mar 2014
#--------------------------------------------------------------------
library(npsp)

# Data
n <- c(15,15)
x1 <- seq(-1,1, length.out = n[1])
x2 <- seq(-1,1, length.out = n[2])

f2d <- function(x,y) x^2 - y^2 
y <- outer(x1, x2, f2d)

# grid.par
range1 <- range(x1)
range2 <- range(x2)
h1 <- diff(range1)/(n[1]-1)
h2 <- diff(range2)/(n[2]-1)

grid <- grid.par(n = n, min = c(min(range1), min(range2)), lag = c(h1, h2))
dim(grid)
str(grid)
coordvalues(grid)
str(coords(grid))  #coords(grid)
dimnames(grid)

# data.grid 
datag <- data.grid(y = y, grid = grid)
dim(datag)
dimnames(datag)
str(datag)
str(coords(datag))  #coords(datag)

spersp(datag, main = 'f(x,y) = x^2 - y^2', xlab = 'x', ylab = 'y')


