# This is file ../spam/tests/xybind.R
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

cat("Testing rbind:\n")
xn <- 3
xm <- 2
yn <- 4
ym <- 2
set.seed(14)

X <- array(runif(xn*xm), c( xn,xm))
Y <- array(runif(yn*ym), c( yn,ym))

R <- as.spam(X)
S <- as.spam(Y)

cat("Testing with two matrices:\n")
test.for.zero( rbind( X, Y), rbind.spam( X, Y))
test.for.zero( rbind( X, Y), rbind.spam( R, Y))
test.for.zero( rbind( X, Y), rbind( R, Y))
test.for.zero( rbind( X, Y), rbind.spam( X, S))
test.for.zero( rbind( X, Y), rbind.spam( R, S))
test.for.zero( rbind( X, Y), rbind( R, S))
test.for.zero( rbind( X, Y*0), rbind.spam( X, as.spam(Y*0)))
test.for.zero( rbind( X*0, Y), rbind.spam( as.spam(X*0), S))
test.for.zero( rbind( X*0, Y*0), rbind.spam( as.spam(X*0), as.spam(Y*0)),rel=F)
test.for.zero( rbind( X*0, Y*0), rbind.spam( X*0, Y*0),rel=F)

cat("Testing with vectors and scalars:\n")
test.for.zero( rbind( X, 1:xm), rbind.spam( X, 1:xm))
test.for.zero( rbind( X, 1:xm), rbind.spam( R, 1:xm))
test.for.zero( rbind( 1:ym, Y), rbind.spam( 1:ym, S))
test.for.zero( rbind( 1, Y), rbind.spam( 1, S))

test.for.zero( rbind( X, 1), rbind.spam( X, 1))

cat("Testing with NULL:\n")
test.for.zero( rbind( X, NULL), rbind.spam( X, NULL))
test.for.zero( rbind( NULL, X, NULL), rbind.spam( NULL, X, NULL))
test.for.zero( rbind( NULL, NULL), rbind.spam( NULL,  NULL),rel=F)

######################################################################

######################################################################
cat("Testing cbind:\n")
xn <- 3
xm <- 2
yn <- 3
ym <- 4
set.seed(14)

X <- array(runif(xn*xm), c( xn,xm))
Y <- array(runif(yn*ym), c( yn,ym))

R <- as.spam(X)
S <- as.spam(Y)

cat("Testing with two matrices:\n")
test.for.zero( cbind( X, Y), cbind.spam( X, Y))
test.for.zero( cbind( X, Y), cbind.spam( R, Y))
test.for.zero( cbind( X, Y), cbind( R, Y))
test.for.zero( cbind( X, Y), cbind.spam( X, S))
test.for.zero( cbind( X, Y), cbind.spam( R, S))
test.for.zero( cbind( X, Y), cbind( R, S))
test.for.zero( cbind( X, Y*0), cbind.spam( X, as.spam(Y*0)))
test.for.zero( cbind( X*0, Y), cbind.spam( as.spam(X*0), S))
test.for.zero( cbind( X*0, Y*0), cbind.spam( as.spam(X*0), as.spam(Y*0)),rel=F)
test.for.zero( cbind( X*0, Y*0), cbind.spam( X*0, Y*0),rel=F)

cat("Testing with vectors and scalars:\n")
test.for.zero( cbind( X, 1:xn), cbind.spam( X, 1:xn))
test.for.zero( cbind( X, 1:xn), cbind.spam( R, 1:xn))
test.for.zero( cbind( 1:yn, Y), cbind.spam( 1:yn, S))
test.for.zero( cbind( 1, Y), cbind.spam( 1, S))

test.for.zero( cbind( X, 1), cbind.spam( X, 1))

cat("Testing with NULL:\n")
test.for.zero( cbind( X, NULL), cbind.spam( X, NULL))
test.for.zero( cbind( NULL, X, NULL), cbind.spam( NULL, X, NULL))
test.for.zero( cbind( NULL, NULL), cbind.spam( NULL,  NULL),rel=F)


######################################################################
options( echo=TRUE)
