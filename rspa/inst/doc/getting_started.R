## ----echo=FALSE, results="hide",include=FALSE----------------------------
library(knitr)
library(rspa)
library(editrules)
opts_chunk$set(size='small')

## ----tidy=FALSE----------------------------------------------------------
E <- editmatrix(expression(
    x5 == x1 + x8,
    x5 == x3 + x4,
    x8 == x6 + x7,
    x4 > 0))

## ----tidy=FALSE----------------------------------------------------------
x <- c(
   x1=330,
   x2=20,
   x3=1000,
   x4=30,
   x5=950,
   x6=500,
   x7=200,
   x8=700)

## ------------------------------------------------------------------------
violatedEdits(E,x,tol=1e-2)

## ------------------------------------------------------------------------
E <- reduce(substValue(E,'x5',x['x5']))

## ------------------------------------------------------------------------
(y <- adjust(E, x))

## ------------------------------------------------------------------------
violatedEdits(E, y$x, tol = 1e-2)

## ----echo=FALSE----------------------------------------------------------
set.seed(1976)

## ----tidy=FALSE----------------------------------------------------------
F <- editmatrix(expression(
   x + y == z,
   x >= 0,
   y >= 0,
   z >=0
))
N <- 100
dat <- data.frame(
   x = runif(100),
   y = rnorm(100),
   z = rlnorm(100)
)

## ------------------------------------------------------------------------
A <- adjustRecords(F,dat)
summary(A)

## ------------------------------------------------------------------------
plot(A)

## ------------------------------------------------------------------------
E

## ----echo=TRUE,tidy=FALSE------------------------------------------------
rc <- data.frame(
   row = c( 1, 1, 2, 2, 3, 3, 3, 4),
   col = c( 1, 2, 3, 4, 2, 5, 6, 4),
  coef = c(-1,-1,-1,-1, 1,-1,-1,-1)
)
b <- c(-950, -950, 0,0)

## ------------------------------------------------------------------------
e <- sparseConstraints(rc, b, neq=3, sorted=TRUE)
e

## ------------------------------------------------------------------------
x_sparse <- c(330, 700, 1000, 30, 500, 200) 
(adjust(e, x_sparse))

## ----eval=FALSE,tidy=FALSE,highlight=FALSE-------------------------------
#  library(LaF)
#  ## Loading required package: Rcpp
#  library(rspa)
#  ## Loading required package: editrules
#  
#  # read A-matrix
#  laf <- laf_open_fwf(
#      file = "prob2A.txt",
#      column_types = c("integer","integer","double"),
#      column_widths = c(10,10,4)
#      )
#  rowcol <- laf[]
#  laf <- close(laf)
#  
#  # read b-vector
#  b <- read.csv("prob2b.txt",header=FALSE)[,1]
#  
#  # read x-vector
#  x <- read.csv("prob2x.txt",header=FALSE)[,1]
#  
#  e <- sparseConstraints(rowcol,b,neq=length(b))
#  e
#  ## Sparse numerical constraints.
#  ##  Variables   : 474948
#  ##  Restrictions: 60675 (printing 0)
#  
#  y <- adjust(e, x)
#  y
#  ## Object of class 'adjusted'
#  ##  Status    : success (using 'sparse' method)
#  ##  Accuracy  : 0.00770226
#  ##  Objective : 47430.5
#  ##  Iterations: 552
#  ##  Timing (s): 5.661
#  ## Solution (truncated at 10):
#  ## [1]  5.4028322  3.5246546  3.7979088  1.1202164  2.5304367  0.2037056
#  ## [7] 97.9161357  1.5714899  4.5743395 -1.1756605

## ------------------------------------------------------------------------
e <- editmatrix(expression( x < 0, x > 1))
isFeasible(e)
adjust(e,c(x=0.5))

