# Reference: 
# Joe H (2006). Generating random correlation matrices based on partial
# correlations. J. Mult. Anal. Vol. 97, 2177--2189.

# Generate random correlation matrix based on random partial correlations
# choice of alpha parameter lead to invariance to index order.
# d = dimension, 
# alphad = alpha parameter for partial of 1,d given 2,...,d-1
# default value alphad = 1 leads to random matrix uniform over
#      space of positive definite correlation matrices
# Other each correlation has a Beta(a,a) distribution on (-1,1) where
#   a=alphad+(d-2)/2
rcorrmatrix<-function(d,alphad=1)
{ 
  d<-as.integer(d)
  if(d<=0 || !is.integer(d))
  { stop("The dimension 'd' should be a positive integer!\n") }
  if(alphad<=0)
  { stop("'alphad' should be positive!\n") }

  # handling of d=1 and d=2
  if(d==1) 
  { rr<-matrix(1,1,1); return(rr) }
  if(d==2) 
  { rho<-runif(1,-1,1)
    rr<-matrix(c(1,rho,rho,1),2,2); return(rr) 
  }
  rr<-matrix(0,d,d)
  diag(rr)<-1
  for(i in 1:(d-1)) 
  { alp<-alphad+(d-2)/2
    rr[i,i+1]<-2*rbeta(1,alp,alp)-1
    rr[i+1,i]<-rr[i,i+1]
  }
  for(m in 2:(d-1))
  { for(j in 1:(d-m))
    { rsub<-rr[j:(j+m),j:(j+m)]
      #print(rsub)
      alp<-alphad+(d-1-m)/2
      rr[j,j+m]<-rjm(rsub,alp)
      rr[j+m,j]<-rr[j,j+m]
    }
  }
  return(rr)
}

# rsub is a symmetrix matrix, alp is a beta parameter
# generate the correlation for the (j,j+m) component
rjm<-function(rsub,alp)
{ b<-nrow(rsub)
  ii<-2:(b-1)
  r1<-rsub[ii,1]
  r3<-rsub[ii,b]
  R2<-rsub[ii,ii]
  Ri<-solve(R2)
  rcond<-2*rbeta(1,alp,alp)-1
  tem13<-t(r1)%*%Ri%*%r3
  tem11<-t(r1)%*%Ri%*%r1
  tem33<-t(r3)%*%Ri%*%r3
  res<-tem13+rcond*sqrt((1-tem11)*(1-tem33))
  return(res)
}



