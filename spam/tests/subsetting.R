# This is file ../spam/tests/subsetting.R
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



# subsetting:
########################################################################


# construct matrices (should be at least 3x5, with n<m):
n <- 10
m <- 15

set.seed(14)
tt <- matrix(rnorm(m*n),n,m)
tt[tt<0] <- 0

ss <- as.spam(tt)

cat("Testing subsetting\n")
test.for.zero(ss[],tt[])      # ok
test.for.zero(ss[,],tt[,])    # ok
test.for.zero(ss[1,],tt[1,])  # ok
test.for.zero(ss[,2],tt[,2])  # ok
test.for.zero(ss[1,3],tt[1,3])# ok
test.for.zero(ss[3:1,],tt[3:1,])# ok


rw <- sample(c(T,F),nrow(tt),rep=T)
cl <- sample(c(T,F),ncol(tt),rep=T)
test.for.zero(ss[rw,cl],tt[rw,cl])
test.for.zero(ss[rw],tt[rw])




rw <- c(1,3);cl <- 1:3;
test.for.zero(ss[rw,cl],tt[rw,cl])
test.for.zero(ss[-rw,cl],tt[-rw,cl])
test.for.zero(ss[-rw,-cl],tt[-rw,-cl])
rw <- c(3,1);cl <- 1:3; test.for.zero(ss[rw,cl],tt[rw,cl])
rw <- c(3,1,2,1);cl <- 1:3; test.for.zero(ss[rw,cl],tt[rw,cl])

tmp <- cbind(sample(1:3,24,rep=T),sample(1:5,24,rep=T))
test.for.zero(ss[tmp],tt[tmp])


test.for.zero(diag(10)[1:2,9:10],diag.spam(10)[1:2,9:10],rel=F)

rs <- sample(c(0,1:(2*n)),2*m,replace=T)
test.for.zero(ss[rs],tt[rs])
# NAs simply work!
rs <- sample(c(0,1:(2*n),NA),2*m,replace=T)
test.for.zero(ss[rs],tt[rs])

rs <- sample(c(T,F,NA),2*m,replace=T)
test.for.zero(ss[rs],tt[rs])

# stuff from 0.31:

tt <- array(1:36,c(6,6))
ss <- as.spam( tt)

for (i in 1:4) {
  rs <- cbind(rep(1:i,each=i),rep(1:i,i))
  test.for.zero(ss[rs],tt[rs])
  test.for.zero(ss[rs+1],tt[rs+1])
  test.for.zero(ss[rs+2],tt[rs+2])

  ti <- array(1:(i^2),c(i,i))
  si <- as.spam( ti)
  test.for.zero(si[rs],ti[rs])  

  si <- spam(0,i,i)
  ti <- as.matrix(si)
  test.for.zero(si[rs],ti[rs])

  si <- diag.spam(i)
  ti <- diag(i)
  test.for.zero(si[rs],ti[rs])
}


if (F) {
# large timing example

  n <- 100
  m <- 15000

  set.seed(14)
  tt <- matrix(rnorm(m*n,mean=-1),n,m)
  tt[tt<0] <- 0

  ss <- as.spam(tt)

  set.seed(14)
  system.time(for (i in 1:100) {
    spam:::subset.rows.spam(ss,sample(1:n,10)) })
  set.seed(14)
  system.time(for (i in 1:100) {
    spam:::subset.spam(ss,sample(1:n,10),1:m) })

}
options( echo=TRUE)
