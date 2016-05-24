el2.test.wtm <- function (xd1,yd1,wxd1new, wyd1new, muvec, nuvec, Hu, Hmu,
  Hnu, p, mean, maxit=10) {

#Initialize vectors and scalars
  lam <- rep(0,p)
  nx1 <- length(muvec)
  ny1 <- length(nuvec)
  swxd1 <- sum(wxd1new)  
  swyd1 <- sum(wyd1new)
  constmat <- matrix(NA, nrow=maxit, ncol=p)
  
for (r in 1:maxit) {
#Calculate muvec1
muvec1 <- rep(NA, nx1)  
  for (k in 1:nx1) {
    muvec1[k] <- wxd1new[k] / abs(swxd1 +  lam %*%
     Hmu[,((k-1)*ny1+1):(k*ny1)] %*% nuvec)
  }
#Calculate nuvec1  
nuvec1 <- rep(NA, ny1)
  for (k in 1:ny1) {
    nuvec1[k] <- wyd1new[k] / abs(swyd1 +  lam %*%
     Hnu[,((k-1)*nx1+1):(k*nx1)] %*% muvec)
  }
#Calculate constraint vector
  constraint <- rep(NA, p)
  for (k in 1:p) {
    constraint[k] <- muvec1 %*% Hu[,((k-1)*ny1+1):(k*ny1)] %*% nuvec1
  }
#Calculate constraintp matrix (derivative of constraint vector wrt lam)
constraintp <- matrix(0, nrow=p, ncol=p)
for (b in 1:p) {
  for (k in 1:p) {
    Hb <- Hu[, ((b-1)*ny1+1):(b*ny1)]
    Hk <- Hu[, ((k-1)*ny1+1):(k*ny1)]
    fact1 <- as.vector((Hk%*%nuvec)*(muvec1^2)/wxd1new)
    fact2 <- as.vector((muvec%*%Hk)*(nuvec1^2)/wyd1new)
    for (i in 1:nx1) {
      for (j in 1:ny1) {
        constraintp[k,b] <- constraintp[k,b] -
        Hb[i,j] * (nuvec1[j]*fact1[i] +
        muvec1[i]*fact2[j]) } } } }
constmat[r,] <- constraint
#Run Newton-Raphson routine
lam1 <- lam - constraint %*% solve(constraintp)
lam <- lam1
}
list(xd1=xd1,wxd1new=wxd1new,muvec1=muvec1,yd1=yd1,wyd1new=wyd1new,
nuvec1=nuvec1,constmat=constmat,lam=lam1,mean=mean)
}
