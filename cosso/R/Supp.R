rescale <- function(x)
  {
  if(length(unique(x))>6)   return( (x-min(x))/(max(x)-min(x)) )
  else                      return(x)
  }

genK<- function(x1,x2)
    {
    m = length(x1)
    n = length(x2)
    b = abs(x1%*%matrix(1,ncol=n)-matrix(1,nrow=m)%*%t(x2))
    k1s = x1 - 0.5
    k1t = x2 - 0.5
    k2s = (k1s^2 - 1/12)/2
    k2t = (k1t^2 - 1/12)/2
    K =k1s%*%t(k1t)+k2s%*%t(k2t)-((b - 0.5)^4 - (b - 0.5)^2 / 2 + 7/240)/24
    }

genK.cat <- function(x1,x2)
  {
  n1 <- length(x1)
  n2 <- length(x2)
  x1 <- rep(x1, times=n2)
  x2 <- rep(x2, each =n1)
  L  <- length(unique(c(x1,x2)))
  K  <- matrix(L*(x1==x2) - 1, n1, n2)
  return(K)
  }

bigGram <- function(x1,x2)
  {
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  d <- ncol(x1)
  Gram <- array(0,c(n1,n2,d))
  for(j in 1:d)
    {
    if(length(unique(c(x1[,j],x2[,j])))>6)   Gram[,,j] <-  genK(    x1[,j], x2[,j])
    else                                     Gram[,,j] <-  genK.cat(x1[,j], x2[,j])
    }
  return(Gram)
  }

wsGram <- function(Gramat,mscale)
   {
    n1 <- dim(Gramat)[1]
    n2 <- dim(Gramat)[2]
    d <- dim(Gramat)[3]
    KK <- matrix(0,n1,n2)
    for (j in 1:d)   KK = KK + mscale[j]*Gramat[,,j]
    return(KK)
   }
   
cvsplitID=function(n,folds)
  {
    fsize  <- floor(n/folds)   #average size of each fold, the last fold is larger
    splits <- fsize*rep(1,folds)
    nextra <- n-folds*fsize
    if (nextra>0)
      {splits[1:nextra] <- splits[1:nextra]+1}
    randid<- sample(1:n,n)
    IDmat <- matrix(NA,ncol=folds,nrow=ceiling(n/folds))
    IDmat[,1]<- randid[1:splits[1]]
    for(i in 2:folds)
       {
        tempid<-randid[(cumsum(splits)[i-1]+1):(cumsum(splits)[i])]
        length(tempid)<- ceiling(n/folds)
        IDmat[,i]<- tempid
       }
  return(IDmat)
  }

My_solve.QP=function(Dmat,dvec,Amat,bvec)
  {
  solution=tryCatch(solve.QP(Dmat,dvec,Amat,bvec)$solution,error=function(x) NA)
  if(is.na(solution[1]))
     {
     Dmat=diag(diag(Dmat))
     solution=solve.QP(Dmat,dvec,Amat,bvec)$solution
     }
   return(solution)
  }
My_solve=function(A,b)
  {
  solution=tryCatch(solve(A,b),error=function(x) NA)
  if(is.na(solution[1]))
     {
     solution=b/diag(A)
     }
   return(solution)
  }

solve.singular=function(A,b)
  {
  solution=tryCatch(solve(A,b),error=function(x) NA)
  if(is.na(solution[1]))
     {
     solution=solve(A+1e-7*diag(nrow(A)),b)
     }
   return(solution)
  }
     