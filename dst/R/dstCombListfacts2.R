dstCombListfacts2<-function(data_table)
  # To combine a serie of support functions with Dempster's rule.
  # Example with three facts
  # a table X of 3 support functions 
  # Note : the complementary mass allocated to Ignorance is not supplied, but added by the function
  # X 
  #    Mass   a    b    c   No  
  #    0.4    0    1    0    1
  #    0.2    1    0    0    1
  #    0.3    1    0    1    3
  #    0.1    0    0    1    3
  #    0.6    1    0    0    5	  
{
  # preeparation of data
  init<-matrix(rep(1,times=ncol(data_table)-1),nrow=1)
  init<-list(DempsterRule=init,con=0)
  data_table<-as.data.frame(data_table)
  n<-nrow(data_table)  
  data_table[,ncol(data_table)]<-as.factor(data_table[,ncol(data_table)])
  nl<-nlevels(data_table[,ncol(data_table)])
  levels(data_table[,ncol(data_table)])<-c(1:nlevels(data_table[,ncol(data_table)]))
  # first pass
  zt<-data_table[,-1]
  zm<-data_table[,c(1,ncol(data_table))]
  card<-ncol(zt)-1
  z1<-subset(zt,zt[,ncol(zt)] == 1)
  if ((nrow(z1) == 1) &  (sum(z1[,-ncol(z1)]) == card)) {
    z1<-init
  } else {
    zm1<-subset(zm,zm[,2] == 1)
    z1a=cbind(zm1[,1],z1[,-ncol(z1)])
    z1b=rbind(z1a,c((1-sum(zm1[,1])), rep(1,times=ncol(z1a)-1)))
    z1<-list(DempsterRule=z1b,con=0)
  }
  zn<-nzdsr(dsrwon(z1,init))
#   Loop
  if (nl == 1) {
    return(zn)
  }  else {
    for (k in 2:nl) {
      zk<-subset(zt,zt[,ncol(zt)] == k)
      if ((nrow(zk) == 1) &  (sum(zk[,-ncol(zk)]) == card)) {
        zk<-init
      } else {
        zm1<-subset(zm,zm[,2] == k)
        z1a=cbind(zm1[,1],zk[,-ncol(zk)])
        z1b=rbind(z1a,c((1-sum(zm1[,1])), rep(1,times=ncol(z1a)-1)))
        zk<-list(DempsterRule=z1b,con=0)      }
      zn<-nzdsr(dsrwon(zn,zk))
    }
    return(zn)
  }
}


