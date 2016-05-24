#This file contains various auxillary functions 

#reorganize constraints matrix according to the predictors and take the difference
#e.g. in a total order the difference is 1
aTx <- function(a, x) {
    if (is.vector(x)) return(x[a[,1]]-x[a[,2]])
    return(x[a[,1],]-x[a[,2],])
}



xT<-function(x) {
   if (is.vector(x)) return(as.matrix(x))
        else return(t(x))
}

#computes A'lambda which should be equal to gradient
taTx<-function(a,x,n) {
  m <- nrow(a)                   
  h <- rep(0,n)                   
  for (i in 1:m) {
    h[a[i,1]] <- h[a[i,1]]+x[i]
    h[a[i,2]] <- h[a[i,2]]-x[i]
  }
return(h)
}

#generates matrix A(I) (i.e. constraint matrix of active sets) as defined in the paper (with 1, -1, and 0)
b2a <- function(b,n) {
  m <- nrow(b)                 #number of active sets
  a <- matrix(0,m,n)
  for (i in 1:m) {             #runs over active sets
    a[i,b[i,1]]<- 1            #1 and -1 depending on constraint specification
    a[i,b[i,2]]<- -1
  }
  return(a)
}

#Warshall algorithm for transitive closure
warshall <- function(a) {
  n <- nrow(a)
  for (j in 1:n) {
    for(i in 1:n) {
      if (a[i,j]==1) a[i,]<-pmax(a[i,],a[j,])
    }
  }
  return(a)
}

#Computes indicators (adjacency, identity, Warshall)
mkIndi<-function(a, n) {
  im <- matrix(0, n, n) 
  m <- nrow(a)
  for (i in 1:m) im[a[i,1],a[i,2]] <- im[a[i,2],a[i,1]] <- 1
  im <- im + diag(n)
  return(t(unique(warshall(im))))
}

#Computes the KKT vector (Lagrange multiplier lambda)
mkLagrange<-function(b, g) {
  ta <- t(b2a(b, length(g)))           #matrix A(I) (1, -1, 0) for active set constraints
  qa <- qr(ta)                         #QR decomposition
  return(qr.coef(qa,g))                #returns KKT vector lambda (solves A*lambda=gradient; stationarity condition)
}

#checks KKT conditions for isotonicity
checkSol<-function(y, gy, a, ay, hl, lbd, ups) {
  ckFeasibility <- min(ay)            #Ax >= 0 (primal feasibility)
  ckLagrange <- min(lbd)              #lambda >= 0 (dual feasibility)
  ckCompSlack <- sum(ay*lbd)          #lambda'Ax = 0 (complementary slackness)
  ckGrad <- max(abs(gy-hl))           #A'lambda-df (stationarity)
return(c(ckGrad,ckFeasibility,ckLagrange,ckCompSlack))
}

weightedMedian<-function(x,w=rep(1,length(x))){
ox<-order(x); x<-x[ox]; w<-w[ox]; k<-1
low<-cumsum(c(0,w)); up<-sum(w)-low; df<-low-up
repeat{
    if (df[k] < 0) k<-k+1
        else if (df[k] == 0) return((w[k]*x[k]+w[k-1]*x[k-1])/(w[k]+w[k-1]))
            else return(x[k-1])
    }
}

weightedFractile<-function(x,w=rep(1,length(x)),a=1,b=1){
ox<-order(x); x<-x[ox]; w<-w[ox]; k<-1
low<-cumsum(c(0,w)); up<-sum(w)-low; df<-a*low-b*up
repeat{
    if (df[k] < 0) k<-k+1
        else if (df[k] == 0) return((w[k]*x[k]+w[k-1]*x[k-1])/(w[k]+w[k-1]))
            else return(x[k-1])
    }
}

weightedMidRange<-function(x,w=rep(1,length(x))){
s<-0; n<-length(x)
if (n==1) return(x)
for (i in 1:(n-1)) for(j in (i+1):n) {
    t<-w[i]*w[j]*abs(x[i]-x[j])/(w[i]+w[j])
    if (t > s) {
        s<-t; i0<-i; j0<-j
        }
    }
return((w[i0]*x[i0]+w[j0]*x[j0])/(w[i0]+w[j0]))
}

#which constraint is active (inactives are 0)
is.active <- function(f,ups=1e-12) which(abs(f) < ups)


is.pos<-function(x,ups=1e-12) x > -ups

is.neg<-function(x,ups=1e-12) x < ups