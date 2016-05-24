# This is file ../spam/tests/crossprod.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


######################################################################
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

spam.options(printsize=60)

######################################################################

cat("Testing crossprod n=1:\n")



set.seed(1)

xf <- rnorm(10)
xf[xf<0] <- 0
xs <- as.spam(xf)

yf <- rnorm(10)
yf[yf<0] <- 0
ys <- as.spam(yf)

test.for.zero( crossprod( xf), crossprod.spam( xs))
test.for.zero( crossprod( xf, yf), crossprod.spam( xs, ys))
test.for.zero( crossprod( xf, yf), crossprod.spam( xs, yf))


dim(xf) <- c(2,5)
dim(yf) <- c(2,5)
ys <- as.spam(yf)
xs <- as.spam(xf)



test.for.zero( crossprod( xf, yf), crossprod.spam( xs, ys))


