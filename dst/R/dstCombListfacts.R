dstCombListfacts<-function(table)
  # To combine a table of facts with Dempster's rule.
  # table : a table of facts. 
  # Example with three facts
  # X  
  #    Mass  oui  non	
  #  	  0.1	  1	  0	  
  #  	  0.1	  0	  1	  
  #  	  0.3	  0	  1	  
{
  init<-matrix(rep(1,times=ncol(table)),nrow=1)
  init<-list(DempsterRule=init,con=0)
  n<-nrow(table)  
  zt<-matrix(table[,2:ncol(table)], nrow=n)
  card<-ncol(zt)
  zm<-table[,1]
  z1<-zt[1,]
  if (sum(z1)==card) {
    z1<-init
  } else {
    z1<-matrix(c(zm[1],z1,(1-zm[1]),rep(1,times=length(z1))),nrow=2, byrow=TRUE)
    z1<-list(DempsterRule=z1,con=0)
  }
  zn<-nzdsr(dsrwon(z1,init))
  if (n == 1) {
    return(zn)
  } else {
    for (k in 2:n) {
      zk<-zt[k,]
      if (sum(zk)==card) {
        zk<-init
      } else {
        zk<-matrix(c(zm[k],zk,(1-zm[k]),rep(1,times=length(zk))),nrow=2, byrow=TRUE)
        zk<-list(DempsterRule=zk,con=0)
      }
      zn<-nzdsr(dsrwon(zn,zk))
    }
    return(zn)
  }
}

