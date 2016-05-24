################################################################################
################################################################################
classif.depth<-function(group,fdataobj,newfdataobj,depth="RP",par.depth=list(),CV="none"){
#,control=list(trace=FALSE,draw=TRUE)
 C<-match.call()
 if (!is.factor(group)) group<-factor(group)
 ismissing<-missing(newfdataobj)
 if (ismissing) newfdataobj<-fdataobj
 group<-factor(group,levels=levels(group)[which(table(group)>0)])
 func.clas<-list()
 lev<-levels(group)
 ng<-length(lev)       
 nc<-ncol(newfdataobj)
 nvec<-table(group)
# p<-nvec[1]/n        
 if (depth %in% c("PD","HD","RP","RPD","RT")){
  if (is.null(par.depth$proj)) {
   d <- nc         
   u <- matrix(runif(d*500,-1,1),500,d)
   norm <- sqrt(rowSums(u*u))
   arg <- u/norm
   if (depth %in% c("RP","RPD","RT"))  par.depth$proj<-fdata(arg,fdataobj$argvals,fdataobj$rangeval)
   else   par.depth$proj<-arg
   } 
 }
 if (depth %in% c("mband","mmode","HD","SD","PD","MhD"))  par.depth$x<-newfdataobj
 else    par.depth[["fdataobj"]]<-newfdataobj 
 depth<-paste("depth.",depth,sep="")

 if (ismissing) {
  ismdist<-is.matrix(par.depth$metric)
 if (ismdist) {
   mdist<-par.depth$metric  
 }
 # if (depth %in% c("HD","SD","PD","MhD"))  par.depth$x<-fdataobj
#  else    par.depth[["fdataobj"]]<-fdataobj 
#  print(names(par.depth))
  n<-nrow(fdataobj)
  x<-array(NA,dim=c(n,nc,ng))
  Df<-matrix(NA,ncol=ng,nrow=n)
#  if (CV!=TRUE){
   ind<-matrix(NA,nrow=n,ncol=ng)
   for (i in 1:ng) {
    ind[,i]<-group==lev[i]
    nam<-c(paste("depth ",lev[i],sep=""),paste("depth ",paste(lev[-i],collapse=",")))
    if (depth %in% c("depth.mband","depth.mmode","depth.SD","depth.HD","depth.PD","depth.MhD"))  par.depth$xx<-fdataobj[ind[,i],]   
    else   par.depth$fdataori<-fdataobj[ind[,i],]  
      if (ismdist) {
          par.depth$metric<-mdist[,ind[,i]]   
          par.depth$metric2<-mdist[ind[,i],ind[,i]]   }     
     Df[,i]<-switch(depth,
     depth.HD= do.call(depth,par.depth)$dep,
     depth.SD= do.call(depth,par.depth)$dep,
     depth.PD= do.call(depth,par.depth)$dep,
     depth.MhD= do.call(depth,par.depth)$dep,
     depth.FM=do.call(depth,par.depth)$dep,                                                                                
     depth.mode=do.call(depth,par.depth)$dep,
     depth.mmode=do.call(depth,par.depth)$dep,    
     depth.RPD=do.call(depth,par.depth)$dep,     
     depth.RP=do.call(depth,par.depth)$dep,
     depth.RT=do.call(depth,par.depth)$dep,
     depth.mband=do.call(depth,par.depth)$dep,      
     depth.band=do.call(depth,par.depth)$dep)
   }
   group.pred<-group.est<-factor(lev[apply(Df,1,which.max)],levels=lev) # Maximum depth 
#  }
  if (CV==TRUE) {
   group.est<-group
   for (j in 1:n) {
    xj<-fdataobj[j,]   
    xnoj<-fdataobj[-j,]   
    ind<-matrix(NA,nrow=n-1,ncol=ng)    
    for (i in 1:ng) {
     ind[,i]<-group[-j]==lev[i]
     xnoji<-xnoj[ind[,i],]    
     nam<-c(paste("depth ",lev[i],sep=""),paste("depth ",paste(lev[-i],collapse=",")))
     if (depth %in% c("depth.mband","depth.mmode","depth.SD","depth.HD","depth.PD","depth.MhD"))  par.depth$xx<-xnoji
     else   par.depth$fdataori<-xnoji
      if (ismdist) {
          par.depth$metric<-mdist[,ind[,i]]   
          par.depth$metric2<-mdist[ind[,i],ind[,i]]   }    
     Df[,i]<-switch(depth,
     depth.HD= do.call(depth,par.depth)$dep,
     depth.PD= do.call(depth,par.depth)$dep,
     depth.SD= do.call(depth,par.depth)$dep,     
     depth.MhD= do.call(depth,par.depth)$dep,
     depth.FM=do.call(depth,par.depth)$dep,
     depth.mode=do.call(depth,par.depth)$dep,
     depth.mmode=do.call(depth,par.depth)$dep,  
     depth.RPD=do.call(depth,par.depth)$dep,     
     depth.RP=do.call(depth,par.depth)$dep,
     depth.band=do.call(depth,par.depth)$dep,
     depth.mband=do.call(depth,par.depth)$dep,      
     depth.RT=do.call(depth,par.depth)$dep)
   } 
  group.est[j]<-factor(lev[which.max(Df[j,])],levels=lev) # Maximum depth }
  }  
  }
  prob.classification<-diag(table(group,group.est))/table(group) 
  mis<-mean(group.est!=group)
  output<-list("group.est"=group.est,"group.pred"=group.pred,"dep"=Df,"depth"=depth, "par.depth"=par.depth,
  "group"=group,"fdataobj"=fdataobj,"C"=C,"prob.classification"=prob.classification,"max.prob"=1-mis)
  class(output)=c("classif")
 return(output) 
}
else  {   # new data
 n<-nrow(newfdataobj)
 n0<-nrow(fdataobj)
 x<-array(NA,dim=c(n,nc,ng))
 Df<-matrix(NA,ncol=ng,nrow=n)
 ind<-matrix(NA,nrow=n0,ncol=ng)
 for (i in 1:ng) {
   ind[,i]<-group==lev[i]
   nam<-c(paste("depth ",lev[i],sep=""),paste("depth ",paste(lev[-i],collapse=",")))
   if (depth %in% c("depth.mband","depth.mmode","depth.SD","depth.HD","depth.PD","depth.MhD"))  par.depth$xx<-fdataobj[ind[,i],]   
   else   par.depth$fdataori<-fdataobj[ind[,i],]  
    Df[,i]<-switch(depth,
    depth.HD= do.call(depth,par.depth)$dep,
    depth.PD= do.call(depth,par.depth)$dep,
    depth.SD= do.call(depth,par.depth)$dep,    
    depth.MhD= do.call(depth,par.depth)$dep,
    depth.FM=do.call(depth,par.depth)$dep,
    depth.mode=do.call(depth,par.depth)$dep,
    depth.mmode=do.call(depth,par.depth)$dep,    
    depth.RP=do.call(depth,par.depth)$dep,
    depth.RPD=do.call(depth,par.depth)$dep,  
    depth.band=do.call(depth,par.depth)$dep,      
    depth.mband=do.call(depth,par.depth)$dep,     
    depth.RT=do.call(depth,par.depth)$dep)
 } 
 group.pred<-factor(lev[apply(Df,1,which.max)],levels=lev) # Maximum depth 
 if (CV!="none"){ 
  if (depth %in% c("mband","mmode","HD","SD","PD","MhD"))  par.depth$x<-fdataobj
  else    par.depth[["fdataobj"]]<-fdataobj 
  n<-nrow(fdataobj)
  x<-array(NA,dim=c(n,nc,ng))
  Df2<-matrix(NA,ncol=ng,nrow=n)
  ind<-matrix(NA,nrow=n,ncol=ng)
  if (!CV) {
   for (i in 1:ng) {
    ind[,i]<-group==lev[i]
    nam<-c(paste("depth ",lev[i],sep=""),paste("depth ",paste(lev[-i],collapse=",")))
    if (depth %in% c("depth.mband","depth.mmode","depth.HD","depth.PD","depth.MhD"))  par.depth$xx<-fdataobj[ind[,i],]   
    else   par.depth$fdataori<-fdataobj[ind[,i],]    
    Df2[,i]<-switch(depth,
    depth.HD= do.call(depth,par.depth)$dep,
    depth.PD= do.call(depth,par.depth)$dep,
    depth.SD= do.call(depth,par.depth)$dep,    
    depth.MhD= do.call(depth,par.depth)$dep,
    depth.FM=do.call(depth,par.depth)$dep,
    depth.mode=do.call(depth,par.depth)$dep,
    depth.mmode=do.call(depth,par.depth)$dep,    
    depth.RP=do.call(depth,par.depth)$dep,
    depth.RPD=do.call(depth,par.depth)$dep,   
    depth.band=do.call(depth,par.depth)$dep,     
    depth.mband=do.call(depth,par.depth)$dep,     
    depth.RT=do.call(depth,par.depth)$dep)
 }
 group.est<-factor(lev[apply(Df2,1,which.max)],levels=lev) # Maximum depth 
 }
 else {
 group.est<-group 
 for (j in 1:n) {                       
   xj<-fdataobj[j,]   
   xnoj<-fdataobj[-j,]  
   ind<-matrix(NA,nrow=n-1,ncol=ng) 
   for (i in 1:ng) {
   ind[,i]<-group[-j]==lev[i]
   xnoji<-xnoj[ind[,i],]    
   nam<-c(paste("depth ",lev[i],sep=""),paste("depth ",paste(lev[-i],collapse=",")))
   if (depth %in% c("depth.mband","depth.mmode","depth.HD","depth.PD","depth.MhD"))  par.depth$xx<-xnoji
   else   par.depth$fdataori<-xnoji

    Df2[,i]<-switch(depth,
    depth.HD= do.call(depth,par.depth)$dep,
    depth.PD= do.call(depth,par.depth)$dep,
    depth.SD= do.call(depth,par.depth)$dep,    
    depth.MhD= do.call(depth,par.depth)$dep,
    depth.FM=do.call(depth,par.depth)$dep,   
    depth.mode=do.call(depth,par.depth)$dep,
    depth.mmode=do.call(depth,par.depth)$dep,    
    depth.RP=do.call(depth,par.depth)$dep,
    depth.RPD=do.call(depth,par.depth)$dep,  
    depth.band=do.call(depth,par.depth)$dep,      
    depth.mband=do.call(depth,par.depth)$dep,          
    depth.RT=do.call(depth,par.depth)$dep)
  } 
  group.est[j]<-factor(lev[which.max(Df2[j,])],levels=lev) # Maximum depth }
  } 
 }
 prob.classification<-diag(table(group,group.est))/table(group) 
 mis<-mean(group.est!=group)
 return(list("group.est"=group.est,"group.pred"=group.pred,"dep"=Df,"dep.ori"=Df2,"depth"=depth,
 "par.depth"=par.depth,"group"=group,"fdataobj"=fdataobj,"C"=C,
 "prob.classification"=prob.classification,"max.prob"=1-mis))
 }
 else return(list("group.pred"=group.pred,"dep"=Df,"depth"=depth,
 "par.depth"=par.depth,"group"=group,"fdataobj"=fdataobj,"C"=C))
   }
 }
################################################################################
################################################################################


