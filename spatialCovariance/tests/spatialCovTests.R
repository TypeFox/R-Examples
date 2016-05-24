library(spatialCovariance)

xVals <- seq(1,10,length=101)
yVals <- NULL

nrows <- 1
ncols <- 1
rowsep <- 0
colsep <- 0

for(x in xVals)
  {
    rowwidth <- x
    colwidth <- x
    info <- precompute(nrows,ncols,rowwidth,colwidth,rowsep,colsep,cat.level=0)
    V <- computeV(info,class="ldt",cat.level=0)
    yVals <- c(yVals,V[1,1])  ## should be 
  }

cat("Passes",sum((yVals-((25-4*pi-4*log(2))/12-log(xVals)))<1e-12),"of 101 tests\n")

nrows <- 2
ncols <- 2
rowsep <- 0
colsep <- 0
rowwidth <- 1
colwidth <- 10
info <- precompute(nrows,ncols,rowwidth,colwidth,rowsep,colsep,cat.level=0)
V1 <- computeV(info,class="ldt",cat.level=0)
rowwidth <- 10
colwidth <- 1
info <- precompute(nrows,ncols,rowwidth,colwidth,rowsep,colsep,cat.level=0)
V2 <- computeV(info,class="ldt",cat.level=0)
cat("Passes", sum(V1 - V2[c(1,3,2,4),c(1,3,2,4)] == 0),"out of 16 additional tests\n")

a <- rowwidth
b <- colwidth
(25*a^2*b^2 - 8 * a* b^3 * atan(a/b) - 8* a^3* b* atan(b/a) - 2*a^4*log(a) -2*b^4*log(b) + (a^4 - 6*a^2*b^2 + b^4)*log(a^2 + b^2))/(12*a^2*b^2)
V2[1,1]