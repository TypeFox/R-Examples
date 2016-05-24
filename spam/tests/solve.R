# This is file ../spam/tests/solve.R
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


# construct spd matrices (should be at least 3x3):
n <- 10


set.seed(11)
tt <- matrix(rnorm(n*n),n,n)
tt <- t(tt) %*% tt
tt[tt<0] <- 0
# I have seen that with R version 2.4.0 Patched (2006-11-25 r39997)
# on i486-pc-linux-gnu, tt is not symmetric...
tt <- tt-(tt-t(tt))/2


ss <- as.spam(tt)

# solving system
cat("Testing 'solve' and derivatives:\n")
b <- rnorm(n)

test.for.zero(solve(ss),solve(tt))
test.for.zero(solve(ss,b),solve(tt,b))



css <- chol(ss)
ctt <- chol(tt[ordering(css),ordering(css)])



test.for.zero(t(as.spam(css))%*%as.spam(css), t(ctt)%*%ctt)
test.for.zero(t(as.spam(css))%*%as.spam(css), tt[ordering(css),ordering(css)])
test.for.zero((t(as.spam(css))%*%as.spam(css))[ordering(css,inv=T),ordering(css,inv=T)], tt)

test.for.zero(backsolve(css,forwardsolve(css,b[ordering(css,inv=T)]))[ordering(css)],
              backsolve(ctt,forwardsolve(t(ctt),b),n))
   #### ,n as patch 

test.for.zero(backsolve(css,b[ordering(css,inv=T)])[ordering(css)],
              backsolve(ctt,b,n))
   #### ,n as patch 

test.for.zero(forwardsolve(css,b[ordering(css,inv=T)])[ordering(css)],
              forwardsolve(t(ctt),b))
test.for.zero(forwardsolve(css,b)[ordering(css)],
              forwardsolve(t(ctt),b[ordering(css)]))

test.for.zero(forwardsolve(css,tt[ordering(css,inv=T),])[ordering(css),],
              forwardsolve(t(ctt),tt))


cat("Testing option 'chol.update' (expect two passes then one fail):\n") 
ss1 <- ss+diag.spam(dim(ss)[1])
test.for.zero( chol(ss), update.spam.chol.NgPeyton(css, ss))

sel <- which(ss[1,,drop=TRUE]!=0)
ss1[1,sel[-1]] <- 0
ss2 <- ss
ss2[n,1] <- .1
spam.options(cholsymmetrycheck=FALSE)
test.for.zero(as.spam(update.spam.chol.NgPeyton(css,ss1)), as.spam( chol(ss1)))
test.for.zero(as.spam(update.spam.chol.NgPeyton(css,ss1)), as.spam( chol(ss2)))
css <- chol(ss)


# spam.options(trivalues=TRUE)
# spam.options(trivalues=FALSE)

spam.options(cholsymmetrycheck=TRUE)

# methods for spam.chol.NgPeyton
cat("Testing methods for 'spam.chol.NgPeyton':\n")
test.for.zero(as.spam(css), ctt)

test.for.zero(as.matrix(css), as.matrix(ctt))
test.for.zero(diag(css), diag(ctt))
test.for.zero(length(css), length(ctt[ctt!=0]))
test.for.zero(dim(css), dim(ctt))
test.for.zero(c(css), c(ctt))


# update singular matrices
cat("Testing option 'cholupdatesingular' (expect a pass, a warning then 2 errors):\n") 
ss3 <- spam(rep(1,4),2)
ch3 <- chol( ss3+diag.spam(2))
spam.options(cholupdatesingular="null")
test.for.zero(is.null(update(ch3, ss3)),TRUE)
spam.options(cholupdatesingular="warning")
options(warn=1)
update(ch3, ss3)
spam.options(cholupdatesingular="error")
try(update(ch3, ss3))
spam.options(cholupdatesingular="NULL")
try(update(ch3, ss3))




# determinants
cat("Testing 'det' and derivatives:\n")
test.for.zero(det(ss),det(tt))
test.for.zero(det(ss,log=T),det(tt,log=T))
test.for.zero(determinant(ss)$mod,determinant(tt)$mod)
test.for.zero(determinant(ss,log=F)$mod,determinant(tt,log=F)$mod)

test.for.zero(det(chol(ss)),det(chol(tt)))

test.for.zero(2*sum(log(diag(css))), determinant(tt)$modulus)


# orderings and derivatives
cat("Testing 'ordering' and derivatives:\n")
tt5 <- matrix(c( 2,0,2,0,4,0,2,0,3),3)
ss5 <- spam(  c( 2,0,2,0,4,0,2,0,3),3)
test.for.zero(ordering(tt5),1:3)
test.for.zero(ordering(ss5),1:3)
test.for.zero(ordering(tt5,inv=T),3:1)
test.for.zero(ordering(ss5,inv=T),3:1)
test.for.zero(ordering(chol(ss5)),c(2,3,1))
test.for.zero(ordering(chol(ss5),inv=T),c(3,1,2))




# spam triangular solves
cat("Testing triangular solves for spam objects:\n")
# We need to generate a upper triangular matrix first.
ctt <- chol(tt)
css <- as.spam(ctt)
b <- rnorm(nrow(tt))


# Recall:
test.for.zero(backsolve(ctt,forwardsolve(t(ctt),b),n),
              solve(tt,b))

# Now do testing:
test.for.zero(forwardsolve(t(css),b), forwardsolve(t(ctt),b))

test.for.zero(forwardsolve(ss,b), forwardsolve(tt,b))

cs <- ss
cs[upper.tri(cs)] <- 0
test.for.zero(forwardsolve(cs,b), forwardsolve(ss,b))


   #### ,n as patch 
test.for.zero(backsolve(css,b), backsolve(ctt,b, n))
test.for.zero(backsolve(ss,b), backsolve(tt,b, n))
test.for.zero(backsolve(t(cs),b), backsolve(tt,b, n))

test.for.zero(backsolve(css,forwardsolve(t(css),b)),
              backsolve(ctt,forwardsolve(t(ctt),b), n))

if (F){
# a few specific tests leading mostly to errors/warnings...

  cs <- css
  cs[3,3] <- 0
  forwardsolve(cs,b)
  backsolve(cs,b)

}

options( echo=TRUE)
