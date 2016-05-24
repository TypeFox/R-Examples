
set.seed(290875)
library("party")

### get rid of the NAMESPACE
attach(asNamespace("party"))

### 
###
###    Regression tests for test statistics
###    
###    functions defined in file `./src/TestStatistic.c'    

### tests for function C_maxabsTeststatistic
xf <- gl(3, 10)
yf <- gl(3, 10)[sample(1:30)]
x <- sapply(levels(xf), function(l) as.numeric(xf == l))
colnames(x) <- NULL
y <- sapply(levels(yf), function(l) as.numeric(yf == l))
colnames(y) <- NULL
weights <- sample(1:30)
linstat <- LinearStatistic(x, y, weights) 
expcov <- ExpectCovarLinearStatistic(x, y, weights)
maxabs <- max(abs(linstat - expcov@expectation) / sqrt(diag(expcov@covariance)))
stopifnot(isequal(maxabs, 
    maxabsTestStatistic(linstat, expcov@expectation, expcov@covariance, 1e-10)))
expcov@covariance[1,1] <- 1e-12
stopifnot(isequal(maxabs,
    maxabsTestStatistic(linstat, expcov@expectation, expcov@covariance, 1e-10)))

### tests for function C_quadformTeststatistic
### -> see LinearStatistic-regtest.R

