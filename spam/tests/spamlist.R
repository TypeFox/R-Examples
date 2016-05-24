# This is file ../spam/tests/spamlist.R
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

test.for.zero(spam( list(ind=numeric(0), j=numeric(0), numeric(0)),nrow=4,ncol=3),
              spam(0,4,3),rel=FALSE)

i <- c(1,2,3,4,5)
j <- c(5,4,3,2,1)
ss3 <- spam(0,5,5)
ss3[cbind(i,j)] <- i/j
test.for.zero(spam.list(list(i=i,j=j,i/j)), ss3)
pad(ss3) <- c(13,13)
test.for.zero(spam.list(list(i=i,j=j,i/j),13,13), ss3)
pad(ss3) <- c(3,3)
test.for.zero(spam.list(list(i=i,j=j,i/j),3,3), ss3)
pad(ss3) <- c(2,2)
test.for.zero(spam.list(list(i=i,j=j,i/j),2,2), ss3,rel=F)


test.for.zero({spam.options(listmethod='EP');
               spam.list(list(i=i,j=j,i/j),ncol=3)},
              {spam.options(listmethod='BS');
               method='BS';spam.list(list(i=i,j=j,i/j),ncol=3)})
test.for.zero({spam.options(listmethod='EP');
               spam.list(list(i=i,j=j,i/j),ncol=3,nrow=4)},
              {spam.options(listmethod='BS');
               spam.list(list(i=i,j=j,i/j),ncol=3,nrow=4)})

test.for.zero(spam.list(list(i=i,j=j,i/j),ncol=1,nrow=1),
              0,rel=F)




set.seed(2011)
m = 1000
rmax = 30
cmax = 40
i = floor(runif(m) * rmax) + 1
j = floor(runif(m) * cmax) + 1
val = floor(10 * runif(m)) + 1


spam.options(listmethod='EP')
ss1 <- spam.list(list(i=i,j=j,val))
spam.options(listmethod='BS')
ss2 <- spam.list(list(i=i,j=j,val))

test.for.zero(ss1,ss2,rel=F)


options( echo=TRUE)
