bonode <-
function(Gmat,y,Xmat,tr,dmatrow,dmatsg,parent,parvec,w,nsplit,crit=crit)  {
   #selects best observed parent node 
     #create matrix to keep the highest value of the criterion
   nleaf<-ncol(Gmat)
   critmax<-matrix(0,nrow=nleaf,ncol=7)
   critmax[,1]<-1:nleaf
   dmatsel<-matrix(0,nrow=dim(dmatsg)[1],ncol=nleaf)
   dmatsel<-sapply(1:nleaf,function(kk,dmatrow,dmatsg,nleaf){ ifelse(apply(dmatsg[,c(nleaf,nleaf+1)]==dmatrow[kk],1,sum)!=2,1,0)},nleaf=nleaf,dmatrow=dmatrow,dmatsg=dmatsg)
   critmax[,c(2:7)]<-t(sapply(1:nleaf,function(kk,y,Xmat,tr,Gmat,dmatsg,dmatsel,parent,parvec,w,nsplit,crit){
     bovar(y,Xmat,tr,Gmat[,kk],dmatsg,dmatsel[,kk],parent[-kk,],parvec,w,nsplit,crit=crit) },
     Xmat=Xmat,Gmat=Gmat,y=y,tr=tr,parent=parent,nsplit=nsplit,w=w,parvec=parvec,dmatsg=dmatsg,dmatsel=dmatsel,crit=crit) )
          #which predictor is the best spllting candidate for this split?
    bestrow<-which(critmax[,5]==max(critmax[,5])) [1]
    #cat("best observed node",critmax[bestrow,],"\n" )
      return(critmax[bestrow,])}
