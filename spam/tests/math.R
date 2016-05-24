# This is file ../spam/tests/math.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


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


# see Matrix::rsparsematrix
spam_random <- function(n, m=n, size=min(m-1,4)*n, fill=rnorm, seed=NULL, ...)
    {
        if (!is.null(seed)) set.seed(seed)
        ind <- sample.int((n*m), size=size)
#
#        as.spam( list(i=(ind %% m)+1, j=(ind %/% n)+1, fill(length(ind), ...)))
        tmp <- matrix(0,n,m)
        tmp[ind] <- fill(length( ind), ...)
        as.spam(tmp)
    }

# construct matrices:
n <- 10
m <- 5

set.seed(14)
tt <- matrix(rnorm(m*n),n,m)
tt[tt<0] <- 0

ss <- as.spam(tt)
spam.options( structurebased=FALSE) # test for equivalence!

#     ‘Math’ ‘"abs"’, ‘"sign"’, ‘"sqrt"’, ‘"ceiling"’, ‘"floor"’,
#          ‘"trunc"’, ‘"cummax"’, ‘"cummin"’, ‘"cumprod"’, ‘"cumsum"’,
#          ‘"log"’, ‘"log10"’, ‘"log2"’, ‘"log1p"’, ‘"acos"’, ‘"acosh"’,
#          ‘"asin"’, ‘"asinh"’, ‘"atan"’, ‘"atanh"’, ‘"exp"’, ‘"expm1"’,
#          ‘"cos"’, ‘"cosh"’, ‘"cospi"’, ‘"sin"’, ‘"sinh"’, ‘"sinpi"’,
#          ‘"tan"’, ‘"tanh"’, ‘"tanpi"’, ‘"gamma"’, ‘"lgamma"’,
#          ‘"digamma"’, ‘"trigamma"’

#     ‘Math2’ ‘"round"’, ‘"signif"’

#     ‘Summary’ ‘"max"’, ‘"min"’, ‘"range"’, ‘"prod"’, ‘"sum"’, ‘"any"’, ‘"all"’

#
# !
A <- diag.spam(4)   ; B <- diag(4)
test.for.zero(A, B)
test.for.zero(!A, !B)
diag(A)=0  ; diag(B) <- 0
test.for.zero(!A, !B)
# str(A) # is what needs to be expected...,
# different to spam:::complement.spam(A)












#     ‘Summary’
test.for.zero(max(ss), max(tt))
test.for.zero(min(ss), min(tt))
test.for.zero(range(ss), range(tt))
test.for.zero(prod(ss), prod(tt))
test.for.zero(sum(ss), sum(tt))
test.for.zero(any(ss), any(tt))
test.for.zero(all(ss), all(tt))

#     ‘Math2’
test.for.zero(round(ss), round(tt))
test.for.zero(signif(ss), signif(tt))

#     ‘Math’ ‘"abs"’, ‘"sign"’, ‘"sqrt"’, ‘"ceiling"’, ‘"floor"’,
#          ‘"trunc"’, ‘"log1p"’
#          ‘"asin"’, ‘"asinh"’, ‘"atan"’, ‘"atanh"’, ‘"expm1"’,
#           ‘"sin"’, ‘"sinh"’, ‘"sinpi"’,
#          ‘"tan"’, ‘"tanh"’, ‘"tanpi"’,

#          ‘"cummax"’, ‘"cummin"’, ‘"cumprod"’, ‘"cumsum"’,
#          ‘"log"’, ‘"log10"’, ‘"log2"’, ‘"acos"’, ‘"acosh"’,
#          , ‘"exp"’, ‘"cos"’, ‘"cosh"’, ‘"cospi"’
#           ‘"gamma"’, ‘"lgamma"’,   ‘"digamma"’, ‘"trigamma"’


test.for.zero(abs(ss), abs(tt))
test.for.zero(cos(ss), cos(tt))
test.for.zero(cosh(ss), cosh(tt))

spam.options( NAOK=TRUE) # test for equivalence!


test.for.zero(gamma(ss), gamma(tt))  #
test.for.zero(digamma(ss), digamma(tt)) #
test.for.zero(trigamma(ss), trigamma(tt))
test.for.zero(exp(ss), exp(tt))
test.for.zero(expm1(ss), expm1(tt))


test.for.zero(log(ss), log(tt))
test.for.zero(cummax(ss), cummax(tt))

for (f in getGroupMembers("Math"))
    test.for.zero( do.call(f, list(ss)),
                  do.call(f, list(tt)), tag=f)
                  




spam.options( structurebased=TRUE) # test for equivalence!
