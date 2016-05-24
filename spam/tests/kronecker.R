# This is file ../spam/tests/kronecker.R
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


spam.options(printsize=6)

xn <- 3
xm <- 2
yn <- 4
ym <- 2
set.seed(14)

X <- array(runif(xn*xm), c( xn,xm))
Y <- array(runif(yn*ym), c( yn,ym))

R <- as.spam(X)
S <- as.spam(Y)

b <- rnorm(5)

# with matrices

cat("Testing with two matrices:\n")
test.for.zero( kronecker( X, Y), kronecker.spam( X, Y) )
test.for.zero( kronecker( X, Y), kronecker.spam( R, S) )


test.for.zero( kronecker( X, Y), kronecker( R, S) )
test.for.zero( kronecker( X, Y), kronecker( R, Y) )
test.for.zero( kronecker( X, Y), kronecker( X, S) )

cat("Testing with a matrix and a vector:\n")
test.for.zero( kronecker( X, b), kronecker.spam( X, b) )
test.for.zero( kronecker( b, Y), kronecker.spam( b, S) )


test.for.zero( kronecker( X, b), kronecker( R, b) )
test.for.zero( kronecker( b, Y), kronecker( b, S) )


cat("Testing degenerate cases\n")
test.for.zero( kronecker( X, 0), kronecker.spam( X, 0),rel=FALSE )
test.for.zero( kronecker( 0, 0), kronecker.spam( 0, 0),rel=FALSE  )
test.for.zero( kronecker( 0, Y), kronecker( spam(0), Y),rel=FALSE  )


cat("Testing for different operators:\n")
test.for.zero( kronecker(X,Y,FUN="+"),kronecker(R,S,FUN="+"))
test.for.zero( kronecker(X,b,FUN="+"),kronecker(R,b,FUN="+"))
test.for.zero( kronecker(c(0,1,0),Y,FUN="+"),kronecker(c(0,1,0),S,FUN="+"))
cat("  expect a warning from testing and from 'kronecker':\n")
test.for.zero( kronecker(diag(2),Y,FUN="+"),kronecker(diag.spam(2),S,FUN="+"))

