# This is file ../spam/tests/ops.R
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

  denom <- ifelse( relative, mean( abs(c(xtrue))),1.0)

  if (any(dim(xtest)!=dim(xtrue)))
    return( cat("## FAILED dimensions  ", dim(xtest), " and ", dim(xtrue),
              fill=TRUE))
  test.value <- sum( abs(c(as.matrix(xtest)) - c( xtrue) ),na.rm=T ) /denom
  if(   test.value < tol ){
          cat("** PASSED test at tolerance ", tol, fill=TRUE)}
  else{ cat( "## FAILED test value = ", test.value, " at tolerance ", tol,
              fill=TRUE)}

}





# construct matrices:
n <- 10
m <- 5

set.seed(14)
tt <- matrix(rnorm(m*n),n,m)
rr <- matrix(rnorm(m*n),n,m)
tt[tt<0] <- 0
rr[rr>0] <- 0

ss <- as.spam(tt)
qq <- as.spam(rr)
spam.options( structurebased=FALSE) # test for equivalence!

spam.options( NAOK=TRUE) # test for equivalence!


for (f in rev(getGroupMembers("Arith")))
    test.for.zero( do.call(f, list(ss,qq)), do.call(f, list(tt,rr)), tag=f)
                  
for (f in getGroupMembers("Compare"))
    test.for.zero( do.call(f, list(ss,qq)), do.call(f, list(tt,rr)), tag=f)
                  
for (f in getGroupMembers("Logic"))
    test.for.zero( do.call(f, list(ss,qq)), do.call(f, list(tt,rr)), tag=f)


tv <- sv <- ss@entries
qv <- qq@entries
spam.options( structurebased=TRUE)


test.for.zero <- function( xtest, xtrue, tol= 1.0e-6, relative=FALSE, tag=NULL){
  # slightly different test function!
  if( !is.null(tag)){
     cat( "testing: ", tag, fill=TRUE)}

  denom <- ifelse( relative, mean( abs(c(xtrue))),1.0)

  test.value <- sum( abs(xtest@entries - c( xtrue) ),na.rm=T ) /denom
  if(   test.value < tol ){
          cat("** PASSED test at tolerance ", tol, fill=TRUE)}
  else{ cat( "## FAILED test value = ", test.value, " at tolerance ", tol,
              fill=TRUE)}

}

for (g in getGroupMembers("Ops")) {
    for (f in getGroupMembers(g)) {
        test.for.zero( do.call(f, list(ss,sv)), do.call(f, list(tv,sv)), tag=f)
        test.for.zero( do.call(f, list(sv,ss)), do.call(f, list(sv,tv)), tag=f)
        test.for.zero( do.call(f, list(ss,4)), do.call(f, list(tv,4)), tag=f)
    }
}
cat("One error caught by 'try':\n")
try(do.call(f, list(ss,1:2)))

####################################################################################################################################

{
spam.options(inefficiencywarning=TRUE)
spam.options(structurebased=FALSE)

diag(2)+diag.spam(2)
}


####################################################################################################################################

options( echo=TRUE)
