# This is file ../spam/tests/overall.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     



options( echo=FALSE)
library( spam, warn.conflict=FALSE)

spam.options(structurebased=FALSE)

test.for.zero <- function( xtest, xtrue, tol= 1.0e-6, relative=TRUE, tag=NULL){

  if( !is.null(tag))     cat( "testing: ", tag, fill=TRUE)

  denom <- ifelse( relative, mean( abs(c(xtrue))),1.0)

  test.value <- sum( abs(c(xtest) - c( xtrue) ) ) /denom
  if( test.value < tol ){  cat("** PASSED test at tolerance ", tol, fill=TRUE)
  }else{ cat( "## FAILED test value = ", test.value, " at tolerance ", tol, fill=TRUE)}

}





# construct matrices (should be at least 3x5, with n<m):
n <- 10
m <- 15

set.seed(14)
tt <- matrix(rnorm(m*n),n,m)
tt[tt<0] <- 0
tt2 <- matrix(rnorm(m*n),m,n)
tt2[tt2<0] <- 0
tt3 <- t(tt2)
tt4 <- tt  %*% t(tt)

ss <- as.spam(tt)
ss2 <- as.spam(tt2)
ss3 <- as.spam(tt3)
ss4 <- ss  %*% t(ss)

########################################################################
# start proper testing...

#ss <- spam(0,5,4)
#ss[cbind(c(1,2,4,5,5),c(1,2,3,2,4))] <- c(1:5)

#4:1 %d+% ss

# add diagonal matrix to spam


########################################################################
kf <- matrix(rnorm(m*n),n,m)
ks <- as.spam(kf)
test.for.zero(max(ks),max(kf))
test.for.zero(range(ks),range(kf))
test.for.zero(log(ks-2*min(ks)),log(kf-2*min(kf)))
test.for.zero(cos(ks),cos(kf))
test.for.zero(round(ks,2),round(kf,2))

########################################################################
cat("Testing transpose\n")
rvec <- 1:n
test.for.zero( t(ss), t(ss))
test.for.zero( t(matrix(rvec,nrow=1)),t(spam(rvec,nrow=1)))
test.for.zero( t(matrix(rvec,ncol=1)),t(spam(rvec,ncol=1)))


########################################################################
# Add/subtract operations:
cat("Testing add/subtracting\n")
test.for.zero(ss+ss,tt+tt)
test.for.zero(ss-ss,tt-tt,rel=FALSE)
test.for.zero(ss+tt,tt+tt)
test.for.zero(tt+ss,tt+tt)
test.for.zero(ss-tt,tt-tt,rel=FALSE)
test.for.zero(tt-ss,tt-tt,rel=FALSE)

test.for.zero(ss+ss3,tt+tt3)
test.for.zero(ss-ss3,tt-tt3)
test.for.zero(ss+tt3,tt+tt3)

test.for.zero(ss+3,tt+3)
test.for.zero(ss+ss[,3,drop=T],tt+tt[,3])
test.for.zero(ss+tt[,3],tt+tt[,3])
test.for.zero(ss+tt[,5],tt+tt[,5])

test.for.zero(3-ss,3-tt)

test.for.zero(ss-tt[,3],tt-tt[,3])
test.for.zero(ss-tt[,5],tt-tt[,5])
test.for.zero(tt[,3]-ss,tt[,3]-tt)
test.for.zero(tt[,5]-ss,tt[,5]-tt)

test.for.zero(-ss,-tt)
test.for.zero(+ss,+tt)

if (F) { # do not run
  ss+t(ss)
  ss+2:(n*m)
  2:(n*m)+ss
  ss-2:(n*m)
  2:(n*m)-ss
  
  ss+ss[,3,drop=F]
  tt+tt[,3,drop=F]
}


# Multiplication
test.for.zero(ss*ss,tt*tt)
test.for.zero(ss*tt,tt*tt)
test.for.zero(tt*ss,tt*tt)


test.for.zero(ss*ss3,tt*tt3)
test.for.zero(ss*ss3,tt*tt3)
test.for.zero(ss*tt3,tt*tt3)
test.for.zero(tt*ss3,tt*tt3)

# Division
test.for.zero(ss*ss,tt*tt)
test.for.zero(ss*tt,tt*tt)
test.for.zero(tt*ss,tt*tt)


test.for.zero(ss*ss3,tt*tt3)
test.for.zero(ss*ss3,tt*tt3)
test.for.zero(ss*tt3,tt*tt3)
test.for.zero(tt*ss3,tt*tt3)


# &,| 
test.for.zero(ss&ss3,tt&tt3)
test.for.zero(ss|ss3,tt|tt3)

# recall, spam has to be first element
test.for.zero(ss&tt3,tt&tt3)
test.for.zero(ss|tt3,tt|tt3)

test.for.zero(ss3&1,tt3&1)
test.for.zero(ss3|1,tt3|1)
test.for.zero(ss3|0,tt3|0)


# subassigning:
cat("Testing assigning\n")
rw <- 1:3
cl <- c(1,3)

ss[1,] <- 1;tt[1,] <- 1
test.for.zero(ss,tt)
ss[1,2] <- 1;tt[1,2] <- 1
test.for.zero(ss,tt)
ss[1,] <- 1:m;tt[1,] <- 1:m
test.for.zero(ss,tt)
ss[3:1,] <- 1:m;tt[3:1,] <- 1:m
test.for.zero(ss,tt)

set.seed(14)
rw <- sample(c(T,F),nrow(tt),rep=T)
cl <- sample(c(T,F),ncol(tt),rep=T)
ass <- rnorm(sum(rw)*sum(cl))
ss[rw,cl] <- ass
tt[rw,cl] <- ass

test.for.zero(ss[rw,cl],tt[rw,cl])

if (F) { # do not run
  # the following quantities will be different because of different reasons.
  # The same is true if we subset only!

  
  # Different row ordering  when assigning
  ss[t(ss2)] <- 1:length(ss2)
  tt[t(tt2)>0] <- 1:length(ss2)

  # Nonzero values are identic but subsetting does not return same structure
  ss[ss2[1:3,1:3]] 
  ss[tt2[1:3,1:3]>0]   # latter two identic
  tt[tt2[1:3,1:3]>0]  



  # The following commands do not work
  ss[ss2] <- 1:length(ss2) # error, wrong format

  ss[as.spam(tmp <- array(1:15,c(5,3)))]  #  error, wrong format

  ss[ array(sample(1:15,24,rep=T),c(12,2))]  # works not because out of bounds

  ss[numeric(0),]  # error, at least one element is needed...
}
ss[cbind(1,1)] <- 4;tt[cbind(1,1)] <- 4
test.for.zero(ss,tt)
ss[rbind(dim(ss))] <- 4;tt[rbind(dim(tt))] <- 4
test.for.zero(ss,tt)

ss[rbind(c(1,1),dim(ss))] <- c(0.1,0.1);tt[rbind(c(1,1),dim(tt))] <- c(0.1,0.1)
test.for.zero(ss,tt)


rw <- c(1,3);cl <- 1:3;
ss[rw,cl] <- 1:3;tt[rw,cl] <- 1:3
test.for.zero(ss,tt)
rw <- c(3,1);cl <- 1:3; 
ss[rw,cl] <- 1:3;tt[rw,cl] <- 1:3
test.for.zero(ss,tt)

nn <- length(tt[-rw,-cl])
test.for.zero(ss[-rw,-cl] <- 1:nn,tt[-rw,-cl] <- 1:nn)
nn <- length(tt[-rw,cl])
test.for.zero(ss[-rw,cl] <- 1:nn,tt[-rw,cl] <- 1:nn)
nn <- length(tt[rw,-cl])
test.for.zero(ss[rw,-cl] <- 1:nn,tt[rw,-cl] <- 1:nn)



# mathstuff
test.for.zero(sqrt(ss),sqrt(tt))


test.for.zero(lower.tri(ss),lower.tri(tt)&tt!=0)
test.for.zero(lower.tri(ss,F),lower.tri(tt,F)&tt!=0)
test.for.zero(upper.tri(ss),upper.tri(tt)&tt!=0)
test.for.zero(upper.tri(ss,F),upper.tri(tt,F)&tt!=0)


if (F) {# only works for full matrices
test.for.zero(ss/tt,tt/tt)
test.for.zero(ss/ss,tt/tt)

kk <- tt/tt
kk[is.na(kk)] <- 0
test.for.zero(ss/tt,kk)
test.for.zero(ss/ss,kk)

test.for.zero(ss^tt,tt^tt)
test.for.zero(ss^ss,tt^tt)
}

# maybe not all of them make sense
if (F) { # this need to be discussed
  ss/ss

  test.for.zero(ss2/tt2,tt2/tt2)
  test.for.zero(ss2^tt2,tt2^tt2)
  test.for.zero(ss2/ss2,tt2/tt2)
  test.for.zero(ss2^ss2,tt2^tt2)
}

# testing rbind/cbind
cat("Testing 'rbind' and 'cbind'\n")
test.for.zero(rbind(tt,t(tt2)), rbind(ss,t(ss2)))
test.for.zero(rbind(tt,tt,t(tt2),1:ncol(tt)), rbind(ss,ss,t(ss2),t(spam(1:ncol(tt)))))

test.for.zero(cbind(tt,t(tt2)), cbind(ss,t(ss2)))
test.for.zero(cbind(tt,tt,t(tt2),1:nrow(tt)), cbind(ss,ss,t(ss2),spam(1:nrow(tt))))

if (F) {
  # dummy testing
  rbind.spam()
  cbind.spam()
  rbind.spam(deparse.level=0)
  cbind.spam(deparse.level=0)
  
  
  # the following should produce warnings:
  rbind(a=ss)
  cbind(b=ss)
  rbind.spam(deparse.level=1)
  cbind.spam(deparse.level=1)



  # the following should produces errors:
  rbind(ss,tt)
  cbind(ss,tt)
  rbind(ss,ss2)
  cbind(ss,ss2)
}


# testing diag:
cat("Testing 'diag' and derivatives:\n")

test.for.zero(diag(tt),diag(ss))
test.for.zero(diag.spam(ss),diag(ss))

test.for.zero(diag.spam(1:4),diag(1:4))
test.for.zero(diag.spam(1,2,3),diag(1,2,3))


test.for.zero(diag.spam(1:4,4,6),diag(1:4,4,6))
test.for.zero(diag.spam(1:4,12),diag(1:4,12))


diag(tt) <- diag(ss) <- 1:n
test.for.zero(tt, ss)
diag(tt) <- diag(ss) <- 2
test.for.zero(tt, ss)

diag(tt) <- diag(ss) <- 0
diag(tt) <- diag(ss) <- 1:n
test.for.zero(tt, ss)

test.for.zero(diag(tt[,5]),diag(ss[,5]))
test.for.zero(diag(tt[,5,drop=T]),diag(ss[,5,drop=T]))
test.for.zero(diag(tt[,5,drop=F]),diag(ss[,5,drop=F]),rel=F)



# testing as.spam
cat("Testing 'as.spam' and derivatives:\n")
test.for.zero(as.spam(0), as.matrix(0),rel=FALSE)
b <- rnorm(n)
test.for.zero(as.spam(b), b )
test.for.zero(as.spam(-abs(b)), -abs(b) )
test.for.zero(-as.spam(abs(b)), -abs(b) )
test.for.zero(as.spam(tt), tt)

# testing spam
test.for.zero(spam(0,1000,1000), matrix(0,1000,1000),rel=FALSE)
test.for.zero(spam(1,12,12),matrix(1,12,12))



# no NA, NaN, Inf handling
if (F) { #not implemented
  ss[2,] <- NA

  tt[4,] <- NA
  tt[,4] <- NaN
  tt[5,] <- Inf
  as.spam(tt)
}



# transformation from list and else
cat("Testing 'as.spam.list' and 'triplet':\n")
test.for.zero( spam(triplet(ss)), ss)
test.for.zero( spam(triplet(ss, tri=TRUE)), ss)

test.for.zero( spam(triplet(tt)), ss)
test.for.zero( spam(triplet(tt, tri=TRUE)), ss)


if (F) { # the following should cause errors:
  spam.list( list(i=1, 2))
  spam.list( list(ind=1, 2))
  spam.list( list(ind=1, 1, 2))
  spam.list( list(ind=1, j=0, 2))
  spam.list( list(ind=numeric(0), j=numeric(0), numeric(0)))
}



options( echo=TRUE)
