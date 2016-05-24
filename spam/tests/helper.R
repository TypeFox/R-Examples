# This is file ../spam/tests/helper.R
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

# bdiag.spam:

A <- spam(rnorm(10),2)
B <- spam(rnorm(16),4)

cat("Testing bdiag.spam:\n")
test.for.zero( bdiag.spam(A),A)

test.for.zero( bdiag.spam(A,B),rbind(cbind(A,rep(0,8)),
                                     cbind(spam(rep(0,20),4),B)))


# rmvnorm:
n <- 5
Sigma <- .25^abs(outer(1:n,1:n,'-'))
Q <- as.spam(solve(Sigma))
b <- 1:n

struct <- chol(Q)
cat("Testing rmvnorm.*:\n")

set.seed(14)
tmp1 <- rmvnorm.canonical(10, b, Q) 
set.seed(14)
test.for.zero( rmvnorm.canonical(10, b, Q, Lstruct=struct), tmp1 )



set.seed(14)
test.for.zero( rmvnorm.prec(10, solve(Q,b), Q), tmp1 )


set.seed(14)
test.for.zero( rmvnorm.prec(10, solve(Q,b), Q, Lstruct=struct), tmp1 )




set.seed(14)
cat("For rmvnorm.canonical:\n- comparing sample mean with truth:\n")
for (i in 10^(1:4))
  cat('    sample size n=',i,' yields  Frobenius-norm:',
      norm( apply(rmvnorm.canonical(i, b, Q, Lstruct=struct), 2,mean)- solve(Q,b),'f'),'\n')
cat("- comparing sample variance with truth:\n")
for (i in 10^(1:4)){
  cat('    sample size n=',i,' yields Frobenius-norm:',
      norm( var( rmvnorm.canonical(i, b, Q=Q, Lstruct=struct))- Sigma,'f'),'\n')
}

set.seed(14)
cat("For rmvnorm.prec:\n- comparing sample mean with truth:\n")
for (i in 10^(1:4))
  cat('    sample size n=',i,' yields  Frobenius-norm:',
      norm( apply(rmvnorm.prec(i, b, Q, Lstruct=struct), 2,mean)- b,'f'),'\n')
cat("- comparing sample variance with truth:\n")
for (i in 10^(1:4)){
  cat('    sample size n=',i,' yields Frobenius-norm:',
      norm( var( rmvnorm.prec(i, Q=Q, Lstruct=struct))- Sigma,'f'),'\n')
}




set.seed(14)
cat("For rmvnorm.spam:\n- comparing sample mean with truth:\n")
for (i in 10^(1:4))
  cat('    sample size n=',i,' yields  Frobenius-norm:',
      norm( apply(rmvnorm.spam(i, b, as.spam(Sigma), Lstruct=struct), 2,mean)- b,'f'),'\n')
cat("- comparing sample variance with truth:\n")
for (i in 10^(1:4)){
  cat('    sample size n=',i,' yields Frobenius-norm:',
      norm( var( rmvnorm.spam(i, b, as.spam(Sigma), Lstruct=struct))- Sigma,'f'),'\n')
}


options( echo=TRUE)
