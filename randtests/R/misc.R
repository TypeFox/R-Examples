#
# miscellaneous functions for R package randtests
#

#
# permut - Generates all permutations of r elements of v
#
# based on function permutations from package gtools
permut <- function (x, m=length(x), FUN=NULL,...){
  n <- length(x)
  X <- NULL  
  if (m == 1) X <- matrix(x, n, 1)
  else if (n == 1) X <- matrix(x, 1, m)
  else if (n == 2 & m == 2) X <- matrix(c(x,x[2:1]), 2, 2)
   else if (n == 4 & m == 4 & is.null(FUN)){
    idx<-c(1:4, c(1:2,4:3), c(1,3,2,4), c(1,3,4,2), c(1,4,2,3), c(1,4,3,2)) 
    X<-rbind(X, matrix(x[idx], nrow=6, ncol=4, byrow=T))
    X<-rbind(X, matrix((x[c(2,1,3:4)])[idx], nrow=6, ncol=4, byrow=T))
    X<-rbind(X, matrix((x[c(3,1,2,4)])[idx], nrow=6, ncol=4, byrow=T))
    X<-rbind(X, matrix((x[c(4,1:3)])[idx], nrow=6, ncol=4, byrow=T))
   }
  else if (n == 5 & m == 5 & is.null(FUN)){
#    X<-matrix(NA,120,5)
    idx<-c(1:5, 1:3,5,4, 1,2,4,3,5, 1,2,4,5,3, 1,2,5,3,4, 1,2,5,4,3, 1,3,2,4,5, 
           1,3,2,5,4, 1,3,4,2,5, 1,3,4,5,2, 1,3,5,2,4, 1,3,5,4,2, 1,4,2,3,5, 
           1,4,2,5,3, 1,4,3,2,5, 1,4,3,5,2, 1,4,5,2,3, 1,4,5,3,2, 1,5,2,3,4, 
           1,5,2,4,3, 1,5,3,2,4, 1,5,3,4,2, 1,5,4,2,3, 1,5,4,3,2)
#    X[1:24,] <- matrix(x[idx], nrow=24, ncol=5, byrow=T)
#    X[25:48,] <- matrix((x[c(2,1,3:5)])[idx], nrow=24, ncol=5, byrow=T)
#    X[49:72,] <- matrix((x[c(3,1,2,4,5)])[idx], nrow=24, ncol=5, byrow=T)
#    X[73:96,] <- matrix((x[c(4,1,2,3,5)])[idx], nrow=24, ncol=5, byrow=T)
#    X[97:120,] <- matrix((x[c(5,1,2,3,4)])[idx], nrow=24, ncol=5, byrow=T)    
     X <- rbind(X, matrix(x[idx], nrow=24, ncol=5, byrow=T))
     X <- rbind(X, matrix((x[c(2,1,3:5)])[idx], nrow=24, ncol=5, byrow=T))
     X <- rbind(X, matrix((x[c(3,1,2,4,5)])[idx], nrow=24, ncol=5, byrow=T))
     X <- rbind(X, matrix((x[c(4,1,2,3,5)])[idx], nrow=24, ncol=5, byrow=T))
     X <- rbind(X, matrix((x[c(5,1,2,3,4)])[idx], nrow=24, ncol=5, byrow=T))
  }  
  else {   
    for (i in 1:n) {
      if(is.null(FUN)) {X <-rbind(X, cbind(x[i], Recall(x[-i], m-1)))}
      else {
        y<-apply(cbind(x[i], Recall(x[-i], m-1)), 1, FUN,...)
        X <-rbind(X,matrix(y))
      }
    }  
  }
  return(X)
} 
#
# original function without the hardcoded cases for 1 < n < 6
#
# permut <- function (x, r=length(x), FUN=NULL,...){ 
#   n <- length(x)
#   X <- NULL  
#   if (r == 1) X<-matrix(x, n, 1)
#   else if (n == 1) X<-matrix(x, 1, r)
#   else {   
#     for (i in 1:n) {
#       if(is.null(FUN)) {X <-rbind(X, cbind(x[i], Recall(x[-i], r-1)))}
#       else {
#         y<-apply(cbind(x[i], Recall(x[-i], r-1)),1,FUN,...)
#         X <-rbind(X,matrix(y))
#       }
#     }  
#   }
#   return(X)
# } 

# randomnesstest - statistic test  
randtests.aux <- function(x, method){
  if (method=="bartels") return(sum(diff(x)^2))
  if (method=="difference.sign"){
    df<-diff(x)
    length(df[df>0])
  } 
}

##
##  probability function of Bartels Rank Statistic
##
dbartelsrank <- function(x, n, log=FALSE){
  stopifnot(is.numeric(x) & n>0)
  tmp <- permut(x=1:n, m=n, FUN=randtests.aux, method="bartels")[,1]
  yr <- rep(0, length(x))
  for (i in 1:length(x)){
    yr[i] <- sum(tmp==x[i])
  }
  r0 <- yr/factorial(n)
  # if TRUE, probabilities p are given as log(p).  
  ifelse(log,return(log(r0)),return(r0))  
}

##
##  distribution function of Bartels Rank Statistic
##
pbartelsrank <- function(q, n, lower.tail=TRUE, log.p=FALSE){
  stopifnot(is.numeric(q) & n>0)
  tmp <- permut(x=1:n, m=n, FUN=randtests.aux, method="bartels")[,1]
  yr <- rep(0, length(q))
  for (i in 1:length(q)){
    yr[i] <- ifelse(lower.tail, sum(tmp<=q[i]), sum(tmp>q[i]))
  }
  r0 <- yr/factorial(n)
  # if TRUE, probabilities p are given as log(p).  
  ifelse(log.p,return(log(r0)),return(r0))
}  