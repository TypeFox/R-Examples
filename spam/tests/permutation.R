# This is file ../spam/tests/permutation.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     








options( echo=FALSE)
library( spam, warn.conflict=FALSE)

test.for.zero <- function( xtest, xtrue, tol= 1.0e-6, relative=TRUE,
tag=NULL){

  if( !is.null(tag)){
     cat( "testing: ", tag, fill=TRUE)}

  denom<-   ifelse( relative, mean( abs(c(xtrue))),1.0)

  test.value <- sum( abs(c(xtest) - c( xtrue) ) ) /denom
  if(   test.value < tol ){
          cat("** PASSED test at tolerance ", tol, fill=TRUE)}
  else{ cat( "## FAILED test value = ", test.value, " at tolerance ", tol,
              fill=TRUE)}

}


########


set.seed(14)
res <- 12.5


grid <- expand.grid(lat=seq(-90+3*res/2,to=90-res,by=res),lon=seq(res/2,to=360,by=res))
dist <- nearest.dist(grid[,2:1],method='gr',upper=NULL, delta=30,R=1)

distm <- as.matrix(dist)

n <- dim(dist)[1]
perm <- sample.int(n,n)

test.for.zero(permutation.spam(dist,P=perm),distm[order(perm),])
test.for.zero(permutation.spam(dist,Q=perm),distm[,order(perm)])
test.for.zero(permutation.spam(dist,P=perm,ind=T),distm[perm,])
test.for.zero(permutation.spam(dist,Q=perm,ind=T),distm[,perm])

test.for.zero(permutation(dist,P=perm),distm[order(perm),])
test.for.zero(permutation(dist,Q=perm),distm[,order(perm)])
test.for.zero(permutation(dist,P=perm,ind=T),distm[perm,])
test.for.zero(permutation(dist,Q=perm,ind=T),distm[,perm])

test.for.zero(permutation(distm,P=perm),distm[order(perm),])
test.for.zero(permutation(distm,Q=perm),distm[,order(perm)])
test.for.zero(permutation(distm,P=perm,ind=T),distm[perm,])
test.for.zero(permutation(distm,Q=perm,ind=T),distm[,perm])

test.for.zero(t(permutation(t(dist),P=perm)),distm[,order(perm)])


