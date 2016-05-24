##### mess by B. Petitpierre, 2011 ####
### mess function calculates the MESS (i.e. extrapolation) as in Maxent ###
### arguments ###
### proj: projection datase t###
### cal: calibration dataset ###
### w: weight for each predictor (e.g. variables importance in SDM) ###
### return values ###
### MESS: mess as calculated in Maxent, i.e. the minimal extrapolation values ###
### MESSw: sum of negative MESS values corrected by the total number of predictors;
### if there is no negative values, MESSw is then the mean MESS ###
### MESSneg: number of predictors on which there is extrapolation ###

ecospat.mess<-function(proj,cal,w="default"){
  if(!is.matrix(proj)){
  proj<-as.matrix(proj)
  cal<-as.matrix(cal)
  }
  if (w=="default"){w<-rep(1,ncol(proj))}

  minp<-apply(cal,2,min);minp<-sapply(minp,rep,t=nrow(proj))
  maxp<-apply(cal,2,max);maxp<-sapply(maxp,rep,t=nrow(proj))
  ecdf.cal<-apply(cal,2,ecdf)

  fn<-function(k,seqdeseq,proj){
  lapply(proj[,k],seqdeseq[[k]])
  }
  fi<-round(unlist(sapply(1:ncol(proj),fn,ecdf.cal,proj,simplify=T))*100)


  #x=fi,proj,minp,maxp
  messi<-function(x){
  if(x[1]==0){MESSi<-(x[2] - x[3]) / (x[4] - x[3]) * 100}
  if(x[1]>0 & x[1]<=50){MESSi<-2*x[1]}
  if(x[1]>50 & x[1]<100){MESSi<-2*(100-x[1])}
  if(x[1]==100){MESSi<-(x[4] - x[2]) / (x[4] - x[3]) * 100}
  return(MESSi)
  }

  count.neg<-function(x){return(length(which(x<0)))}
  total<-round(matrix(apply(cbind(fi,as.vector(proj),as.vector(minp),as.vector(maxp)),1,messi),nrow=nrow(proj)))
  if(ncol(proj)>1){
  MESS<-apply(total,1,min)
  MESSneg<-total;MESSneg[MESSneg[]>0]<-0
  MESSneg[which(MESS>=0),]<-total[which(MESS>=0),]
  MESSw<-round(apply(MESSneg,1,weighted.mean,w=w))
  MESSneg<-apply(total,1,count.neg)
  return(cbind(MESS,MESSw,MESSneg))}else{
  return(total)
  }
}

##### plot.mess by B. Petitpierre, 2011 ####
### plot the MESS extrapolation index onto the geographical space ###
### arguments ###
### xy: xy coordinates of the projection dataset t###
### mess.object: dataframe returned by the mess() function ###
### return values ###
### MESS: mess as calculated in Maxent, i.e. the minimal extrapolation values
### (red= negative, blue= positive values) ###
### MESSw: sum of negative MESS values corrected by the total number of predictors;
### if there is no negative values, MESSw is then the mean MESS (red= negative,
### blue= positive values)###
### MESSneg: number of predictors on which there is extrapolation ###

ecospat.plot.mess<-function(xy,mess.object,cex=1,pch=15){

col.mess.neg<-colorRampPalette(c("white","red"))
col.mess.pos<-colorRampPalette(c("white","blue"))

col.neg<-col.mess.neg(max(1+abs(mess.object[,1])))
col.pos<-col.mess.pos(max(1+abs(mess.object[,1])))
par(mfrow=c(2,2))
#x11()
plot(xy,cex=cex,pch=pch,col=0, main="MESS",xlab=paste("min=",min(mess.object[,1])," & max=",max(mess.object[,1]),sep=""),ylab="")
points(xy[which(mess.object[,1]<0),],cex=cex,pch=pch,col=col.neg[mess.object[which(mess.object[,1]<0),1]])
points(xy[which(mess.object[,1]>0),],cex=cex,pch=pch,col=col.pos[mess.object[which(mess.object[,1]>0),1]])

col.neg<-col.mess.neg(max(1+abs(mess.object[,2])))
col.pos<-col.mess.pos(max(1+abs(mess.object[,2])))
#x11()
plot(xy,cex=cex,pch=pch,col=0, main="MESSw",xlab=paste("min=",min(mess.object[,2])," & max=",max(mess.object[,2]),sep=""),ylab="")
points(xy[which(mess.object[,2]<0),],cex=cex,pch=pch,col=col.neg[mess.object[which(mess.object[,2]<0),2]])
points(xy[which(mess.object[,2]>0),],cex=cex,pch=pch,col=col.pos[mess.object[which(mess.object[,2]>0),2]])

#x11()
col.neg<-col.mess.neg(max(1+abs(mess.object[,3])))
plot(xy,cex=cex,pch=pch,col=col.neg[mess.object[,3]+1], main="#MESSneg",xlab=paste("min=",min(mess.object[,3])," & max=",max(mess.object[,3]),sep=""),ylab="")
}
