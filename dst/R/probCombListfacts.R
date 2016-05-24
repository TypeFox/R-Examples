probCombListfacts<-function(table)
  # to combine a table of facts on a {0, 1} (yes, no) frame, without ignorance
  # (probability distributions) with Dempster's rule  
  # Example of a table
  # table
  #       [,1] [,2] [,3] 
  # [1,]  0.1    1    0
  # [2,]  0.1    0    1
  # [3,]  0.3    0    1
{
 # init
  init<-list(DempsterRule=matrix(c(1,1,1),nrow=1),con=0)
  zt<-as.matrix(table[,2:3])
  zm<-table[,1]
  n<-nrow(zt)
  z1<-zt[1,]
  z1<-matrix(c(zm[1],z1,(1-zm[1]),(c(1,1)-z1),0,1,1),nrow=3, byrow=TRUE)
  z1<-list(DempsterRule=z1,con=0)
  zn<-nzdsr(dsrwon(z1,init))
  # loop
  for (k in 2:n)
  {
    zk<-zt[k,]
    zk<-matrix(c(zm[k],zk,(1-zm[k]),(c(1,1)-zk),0,1,1),nrow=3, byrow=TRUE)
    zk<-list(DempsterRule=zk,con=0)
    zn<-nzdsr(dsrwon(zn,zk))
  }
  return(zn)
}

