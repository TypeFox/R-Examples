## ------------------------------------------------------------------------
library(ddR)
useBackend(parallel,executors=2)
## Run the next two lines also to use the Distributed R backend
# library(distributedR.ddR)
# useBackend(distributedR)

## ------------------------------------------------------------------------
my.dlist <- dlist(1,2,3,4,5)
my.dlist

## ------------------------------------------------------------------------
my.dlist <- dlist(1,2,3,4,5,nparts=3)
my.dlist

## ------------------------------------------------------------------------
my.darray <- darray(dim=c(4,4),psize=c(2,2),data=3)
my.darray

## ------------------------------------------------------------------------
my.dlist <- dlapply(1:5,function(x) x)
my.dlist

## ------------------------------------------------------------------------
my.dlist <- dlapply(1:5,function(x) x, nparts=3)
my.dlist

## ------------------------------------------------------------------------
my.darray2 <- dmapply(function(x) matrix(x,2,2), 1:4, output.type="darray", combine="rbind", nparts=c(2,2))
my.darray2

## ------------------------------------------------------------------------
my.array <- collect(my.darray2)
my.array

## ------------------------------------------------------------------------
collect(my.darray2,3)

## ------------------------------------------------------------------------
my.darray2

## ------------------------------------------------------------------------
parts(my.darray2)

## ------------------------------------------------------------------------
parts(my.darray2,2:3)

## ----eval=FALSE----------------------------------------------------------
#  dlist1 <- dlapply(arg,FUN,nparts)
#  dlist2 <- dmapply(FUN,arg1,arg2,MoreArgs,nparts)
#  darray.or.dframe <- dmapply(FUN,arg1,arg2,MoreArgs,output.type,combine,nparts)

## ------------------------------------------------------------------------
## Head and tail
head(my.darray2,n=1)
tail(my.darray2,n=1)

## Subsetting
my.darray2[2,c(2,1)]

## Statistics
colSums(my.darray2)
max(my.darray2)

## ------------------------------------------------------------------------
da <- darray(psize=c(2,2),dim=c(4,4),data=3)
da

## ------------------------------------------------------------------------
skel <- darray(psize=c(4,2),dim=c(4,4),data=0)

## ------------------------------------------------------------------------
da <- repartition(da,skel)
da

## ------------------------------------------------------------------------
collect(da)

## ------------------------------------------------------------------------
a <- dmapply(function(x) { x }, rep(3,5))
collect(a)

## ------------------------------------------------------------------------
a

## ------------------------------------------------------------------------
b <- dmapply(function(x,y) { x + y }, a, 1:5,nparts=1)
b

## ------------------------------------------------------------------------
collect(b)

## ------------------------------------------------------------------------
addThenSubtract <- function(x,y,z) {
  x + y - z
}
c <- dmapply(addThenSubtract,a,b,MoreArgs=list(z=5))
collect(c)

## ------------------------------------------------------------------------
d <- dmapply(function(x) length(x),parts(a))
collect(d)

## ------------------------------------------------------------------------
e <- dmapply(function(x) length(x),parts(b))
collect(e)

## ----eval=FALSE----------------------------------------------------------
#  library(distributedR.ddR)

## ----eval=FALSE----------------------------------------------------------
#  useBackend(distributedR,executors=2)

