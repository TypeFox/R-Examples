## ------------------------------------------------------------------------
  matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2, byrow = TRUE)

## ------------------------------------------------------------------------
matrix(c(1, 2, 3, 4,
         5, 6, 7, 8), nrow = 2, byrow = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("ramify")  # latest stable release

## ------------------------------------------------------------------------
library(ramify)  # load package
mat("1, 2, 3, 4; 5, 6, 7, 8")

## ------------------------------------------------------------------------
mat("1:4; 5:8", rows = FALSE)  # ; separates columns

## ------------------------------------------------------------------------
mat("1 2 3 4; 5 6 7 8", sep = "")  # blank spaces separate columns

## ------------------------------------------------------------------------
options(mat.rows = FALSE, mat.sep = "")  # change default behavior

## ------------------------------------------------------------------------
strsplit("rnorm(10, sd = 3)", split = ",")

## ------------------------------------------------------------------------
mat("rnorm(5); rnorm(5, mean = 10)", sep = NULL)

## ------------------------------------------------------------------------
z <- list("a" = 1:10, "b" = 11:20, "c" = 21:30)

## ------------------------------------------------------------------------
# Approach 1
matrix(unlist(z), nrow = 3, byrow = TRUE)

# Approach 2
do.call(rbind, z)

# Approach 3
t(simplify2array(z))

## ------------------------------------------------------------------------
mat(z) 

## ------------------------------------------------------------------------
m <- matrix(rnorm(1000000), nrow = 1000)
# head(m)

## ------------------------------------------------------------------------
m <- randn(1000, 1000)  # see Table 2 for a description of randn
pprint(m)

## ------------------------------------------------------------------------
# List holding individual variables
z1 <- list(Height = c(Joe = 6.2,   Mary = 5.7,   Pete = 6.1),
           Weight = c(Joe = 192.2, Mary = 164.3, Pete = 201.7),
           Gender = c(Joe = 0,     Mary = 1,     Pete = 0))

as.data.frame(z1)  # convert z1 to a data frame

## ------------------------------------------------------------------------
# List holding records (i.e., individual observations)
z2 <- list(Joe  = c(Height = 6.2, Weight = 192.2, Gender = 0),
           Mary = c(Height = 5.7, Weight = 164.3, Gender = 1),
           Pete = c(Height = 6.1, Weight = 201.7, Gender = 0))

## ------------------------------------------------------------------------
as.data.frame(t(as.data.frame(z2)))

## ------------------------------------------------------------------------
dmat(z1, rows = FALSE)  # treat list elements as columns of a data frame
dmat(z2)  # treat list elements as rows of a data frame

## ------------------------------------------------------------------------
A1 <- matrix(c(1, 2, 5, 6), nrow = 2, byrow = TRUE)
A2 <- matrix(c(3, 4, 7, 8), nrow = 2, byrow = TRUE)
A3 <- matrix(c(9, 10, 11, 12), nrow = 1)
A <- rbind(cbind(A1, A2), A3)

## ------------------------------------------------------------------------
# A1 <- mat("1, 2; 5, 6")
# A2 <- mat("3, 4; 7, 8")
# A3 <- mat("9, 10, 11, 12")
# A <- bmat("A1, A2; A3")

## ------------------------------------------------------------------------
m <- rand(10, 10)
sum(diag(m))  # compute the trace in base R
tr(m)  # compute the trace using ramify's tr function

## ------------------------------------------------------------------------
a <- array(runif(20000), c(100, 100, 2))  # base R
a <- rand(100, 100, 2)  # ramify
pprint(a[,,1])  # print the first matrix

## ------------------------------------------------------------------------
zeros(10)            # 10x1 matrix (i.e., column vector) of zeros
ones(10, 10)         # 10x10 matrix of ones
# fill(pi, 10, 10, 3)  # 10x10x3 array filled with the constant pi

## ------------------------------------------------------------------------
zeros(10, atleast_2d = FALSE)  # has class "integer", rather than "matrix"

## ---- fig.cap='Color image and contour lines for $y$.', fig.width=7, fig.height=5----
x <- meshgrid(linspace(-4*pi, 4*pi, 27))  # list of input matrices
y <- cos(x[[1]]^2 + x[[2]]^2) * exp(-sqrt(x[[1]]^2 + x[[2]]^2)/6)
par(mar = c(0, 0, 0, 0))  # remove margins
image(y, axes = FALSE)  # color image
contour(y, add = TRUE, drawlabels = FALSE)  # add contour lines

