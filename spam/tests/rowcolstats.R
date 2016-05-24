# This is file ../spam/tests/rowcolstats.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     








#options( echo=FALSE)
library( spam, warn.conflict=FALSE)


test.for.zero <- function( xtest, xtrue, tol= 1.0e-6, relative=FALSE,
tag=NULL){
  # slightly different test function!
  if( !is.null(tag)){
     cat( "testing: ", tag, fill=TRUE)}

  denom<-   ifelse( relative, mean( abs(c(xtrue))),1.0)

  if (any(dim(xtest)!=dim(xtrue)))
    return( cat("## FAILED dimensions  ", dim(xtest), " and ", dim(xtrue),
              fill=TRUE))
  test.value <- sum( abs(c(xtest) - c( xtrue) ),na.rm=T ) /denom
  if(   test.value < tol ){
          cat("** PASSED test at tolerance ", tol, fill=TRUE)}
  else{ cat( "## FAILED test value = ", test.value, " at tolerance ", tol,
              fill=TRUE)}

}



# simple tests:
########################################################################


# construct matrices:
n <- 10
m <- 15

set.seed(14)
tt <- matrix(rnorm(m*n),n,m)
tt[tt<0] <- 0

ss <- as.spam(tt)

test.for.zero(rowSums(ss),rowSums(tt)) 
test.for.zero(colSums(ss),colSums(tt)) 


spam.options(structurebased=FALSE)

test.for.zero(rowMeans(ss),rowMeans(tt))      # ok
test.for.zero(colMeans(ss),colMeans(tt))      # ok



spam.options(structurebased=TRUE)
test.for.zero(rowMeans(ss),rowSums(tt)/apply(tt>0,1,sum))      # ok
test.for.zero(colMeans(ss),colSums(tt)/apply(tt>0,2,sum))      # ok

test.for.zero(rowMeans(ss),apply.spam(ss,1,mean))      # ok
test.for.zero(colMeans(ss),apply.spam(ss,2,mean))      # ok


test.for.zero(rowMeans(spam(0,n,m)),rowMeans(tt*0))      # ok
test.for.zero(colMeans(spam(0,n,m)),colMeans(tt*0))      # ok
test.for.zero(rowMeans(as.spam(diag(0,n))),rowMeans(diag(0,n)))      # ok
test.for.zero(colMeans(as.spam(diag(0,n))),colMeans(diag(0,n)))      # ok

spam.options(structurebased=TRUE)


options( echo=TRUE)
