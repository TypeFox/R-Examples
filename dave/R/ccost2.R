ccost2 <-
function(veg,oldgr,newgr,y=1) {
# call:  o.ccost<- ccost(nveg,oldgr,newgr,y=0.25)
# computing cost function from data frame veg,
# old and new classification with transformation
# coefficient y. 
 nrel<- length(veg[,1])
 ng<- length(table(oldgr))
 C<- matrix(c(rep(0,ng*ng)),ncol=ng)
 for(i in 1:nrel) C[oldgr[i],newgr[i]]<- C[oldgr[i],newgr[i]]+1
 C<- as.matrix(C,ncol=ng)
 C<- as.matrix(C,ncol=ng)
# rownames(C) <- table(oldgr)
# colnames(C) <- table(oldgr)
 order.r<- seq(1,nrel,1)
 o.cr<- rcentroid(veg,oldgr,y)
 W<- as.matrix(o.cr$dist.mat,ncol=ng)
 cf<- sum(C*W)
# cat("cf=",cf,"\n")
 o.ccost<- list(dimension=ng,ccost=cf,old.groups=oldgr,new.groups=newgr,conf.matrix=C,weight.matrix=W,transf=y)
 }
