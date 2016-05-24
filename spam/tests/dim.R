# This is file ../spam/tests/dim.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


options( echo=FALSE)
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
  test.value <- sum( abs(c(as.matrix(xtest)) - c( xtrue) ),na.rm=T ) /denom
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




test.for.zero(ss,tt)

dim(ss) <- c(m,n)
dim(tt) <- c(m,n)
test.for.zero(ss,tt)

dim(ss) <- c(m*n,1)
dim(tt) <- c(m*n,1)
test.for.zero(ss,tt)

dim(ss) <- c(1, m*n)
dim(tt) <- c(1, m*n)
test.for.zero(ss,tt)

cat("Two obvious errors caught by 'try':\n")
try( dim(ss) <- c(-1, -m*n))
try( dim(ss) <- c(1, m, n))


options( echo=TRUE)
