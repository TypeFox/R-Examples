rmanova <- function(y, groups, blocks, tr = 0.2){
  
  cols1 <- deparse(substitute(y))
  cols2 <- deparse(substitute(groups))
  cols3 <- deparse(substitute(blocks))
  dat <- data.frame(y, groups, blocks)
  colnames(dat) <- c(cols1, cols2, cols3)
  cl <- match.call()
 
  x <- reshape(dat, idvar = cols3, timevar = cols2, direction = "wide")[-1]  ## wide format
  grp <- c(1:length(x))
  
  if(is.list(x)){
    J<-length(grp)  # The number of groups to be compared
    m1<-matrix(x[[grp[1]]],length(x[[grp[1]]]),1)
    for(i in 2:J){     # Put the data into an n by J matrix
      m2<-matrix(x[[grp[i]]],length(x[[i]]),1)
      m1<-cbind(m1,m2)
    }
  }
  
  m2<-matrix(0,nrow(m1),ncol(m1))
  xvec<-1
  g<-floor(tr*nrow(m1))  #2g is the number of observations trimmed.
  for(j in 1:ncol(m1)){  # Putting Winsorized values in m2
    m2[,j]<-winval(m1[,j],tr)
    xvec[j]<-mean(m1[,j],tr)
  }
  xbar<-mean(xvec)
  qc<-(nrow(m1)-2*g)*sum((xvec-xbar)^2)
  m3<-matrix(0,nrow(m1),ncol(m1))
  m3<-sweep(m2,1,apply(m2,1,mean))  # Sweep out rows
  m3<-sweep(m3,2,apply(m2,2,mean))  # Sweep out columns
  m3<-m3+mean(m2)  # Grand Winsorized mean swept in
  qe<-sum(m3^2)
  test<-(qc/(qe/(nrow(m1)-2*g-1)))
  #
  #  Next, estimate the adjusted degrees of freedom
  #
  v<-winall(m1,tr=tr)$cov
  vbar<-mean(v)
  vbard<-mean(diag(v))
  vbarj<-1
  for(j in 1:J){
    vbarj[j]<-mean(v[j,])
  }
  A<-J*J*(vbard-vbar)^2/(J-1)
  B<-sum(v*v)-2*J*sum(vbarj^2)+J*J*vbar^2
  ehat<-A/B
  etil<-(nrow(m2)*(J-1)*ehat-2)/((J-1)*(nrow(m2)-1-(J-1)*ehat))
  etil<-min(1.,etil)
  df1<-(J-1)*etil
  df2<-(J-1)*etil*(nrow(m2)-2*g-1)
  siglevel<-1-pf(test,df1,df2)
  result <- list(test=test,df1 = df1, df2 = df2, p.value = siglevel, call = cl)
  class(result) <- c("t1way")
  result
}
