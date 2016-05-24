dlproj<-function(dlobj=NULL,projyears=NULL,
  projtype=1,projcatch=NULL,grout=1,
  grargs = list(lwd=1,unit="MT",mains=" ",cex.main=1,cex.axis=1,cex.lab=1),
  grtif=list(zoom=4,width=11,height=13,pointsize=10))
  {
  if(is.null(dlobj)) stop("dlobj is missing")
  if(!any(names(dlobj)=="type")) stop("dlobj not a catchmsy or dbsra output object file")
  if(!dlobj$type %in% c("dbsra","catchmsy")) stop("dlobj not a catchmsy or dbsra output object file")
  if(dlobj$type %in% c("dbsra","catchmsy")){
  
    if(projtype==2 & is.null(projcatch)) stop("For projtype=2, enter value for projcatch")
    if(is.null(projyears)) stop("projyears is missing")
  }
  
  #gr def
  grdef=list(lwd=1,unit="MT",cex.axis=1,cex.lab=1,
        mains=" ",cex.main=1)
   if(any(names(grargs)=="unit")){ 
            if(any(grdef$unit!=grargs$unit)) grdef$unit<-grargs$unit
   }
  if(any(names(grargs)=="cex.axis")){ 
            if(any(grdef$cex.axis!=grargs$cex.axis)) grdef$cex.axis<-grargs$cex.axis
   }
   if(any(names(grargs)=="cex.lab")){ 
            if(any(grdef$cex.lab!=grargs$cex.lab)) grdef$cex.lab<-grargs$cex.lab
   }
   if(any(names(grargs)=="mains")){ 
            if(any(grdef$mains!=grargs$mains)) grdef$mains<-grargs$mains
   }
   if(any(names(grargs)=="cex.main")){ 
            if(any(grdef$cex.main!=grargs$cex.main)) grdef$cex.main<-grargs$cex.main
   }
   if(any(names(grargs)=="lwd")){ 
            if(any(grdef$lwd!=grargs$lwd)) grdef$lwd<-grargs$lwd
   }
   tifdef=list(zoom=4,width=11,height=13,pointsize=10)
    if(any(names(grtif)=="zoom")){ 
            if(tifdef$zoom!=grtif$zoom) tifdef$zoom<-grtif$zoom
   }
  if(any(names(grtif)=="width")){ 
            if(tifdef$width!=grtif$width) tifdef$width<-grtif$width
   }
  if(any(names(grtif)=="height")){ 
            if(tifdef$height!=grtif$height) tifdef$height<-grtif$height
   }       
   if(any(names(grtif)=="pointsize")){ 
            if(tifdef$pointsize!=grtif$pointsize) tifdef$pointsize<-grtif$pointsize
   }   
  
################ dbsra ##################################
  if(dlobj$type=="dbsra"){
      agemat<-dlobj$agemat
     #projtype(0=median MSY, 1=mean MSY, 2=User input) 
      if(projtype<2){ 
       if(projtype==0) catch1<-dlobj$Estimates[1,2] #use median as catch 
       if(projtype==1) catch1<-dlobj$Estimates[1,1] #use mean
       catch1<-rep(catch1,projyears)
      } 
   if(projtype==2){# user input
      if(is.null(projcatch)) stop("projcatch must be specified")
     if(length(projcatch)>1 && length(projcatch)<projyears) stop("length of projcatch does not match length of projyears")
     if(length(projcatch)>1)  catch1<-projcatch
     if(length(projcatch)==1) catch1<-rep(projcatch,projyears)
   }
    
    parms<-dlobj$Values[dlobj$Values[,1]==1,c(2,4,5,6)]
    names(parms)<-c("FmsyM","BmsyK","M","K") 
     w<-function(x){
        f<-function(d){(x-(d^(1/(1-d))))^2}
        optimize(f,c(0,1000),tol=0.000001)[[1]]
     }
     parms$n<-apply(as.data.frame(parms$BmsyK),1,w)
     parms$g<-(parms$n^(parms$n/(parms$n-1)))/(parms$n-1)
     P<-NULL
     bigB<-read.csv("Biotraj-dbsra.csv",header=FALSE)[,-1]
     bigB<-t(bigB)    
       
    pB<-as.data.frame(bigB[c(nrow(bigB)-agemat):nrow(bigB),c(which(dlobj$Values[,1]==1))])
    star<-nrow(pB)
    pB[c(nrow(pB)+1):c(c(nrow(pB)+1)+projyears-1),]<-NA
    
     for(nn in 1:nrow(parms)){
        P<-0
        Fmsy<-parms$FmsyM[nn]*parms$M[nn]
        Umsy<-(Fmsy/(Fmsy+parms$M[nn]))*(1-exp(-Fmsy-parms$M[nn]))
        MSY<-parms$K[nn]*parms$BmsyK[nn]*Umsy
      
       for(t in star:c(nrow(pB)-1)){  
           if(parms$BmsyK[nn]>=0.5) P<-parms$g[nn]*MSY*(pB[t-agemat,nn]/parms$K[nn])-parms$g[nn]*MSY*(pB[t-agemat,nn]/parms$K[nn])^parms$n[nn] 
           if(parms$BmsyK[nn]>0.3 & parms$BmsyK[nn]<0.5){
              bjoin<-(0.75*parms$BmsyK[nn]-0.075)*parms$K[nn]
             if(pB[t-agemat,nn]<bjoin){
               PJ<-parms$g[nn]*MSY*(bjoin/parms$K[nn])-parms$g[nn]*MSY*(bjoin/parms$K[nn])^parms$n[nn]
               cc<-(1-parms$n[nn])*MSY*parms$g[nn]*(bjoin^(parms$n[nn]-2))*parms$K[nn]^-parms$n[nn]
               P<-pB[t-agemat,nn]*(PJ/bjoin+cc*(pB[t-agemat,nn]-bjoin)) 
              }
             if(pB[t-agemat,nn]>=bjoin) P<-parms$g[nn]*MSY*(pB[t-agemat,nn]/parms$K[nn])-parms$g[nn]*MSY*(pB[t-agemat,nn]/parms$K[nn])^parms$n[nn] 
           }
           if(parms$BmsyK[nn]<=0.3){
              bjoin<-(0.5*parms$BmsyK[nn])*parms$K[nn]
            if(pB[t-agemat,nn]<bjoin){
               PJ<-parms$g[nn]*MSY*(bjoin/parms$K[nn])-parms$g[nn]*MSY*(bjoin/parms$K[nn])^parms$n[nn]
               cc<-(1-parms$n[nn])*MSY*parms$g[nn]*(bjoin^(parms$n[nn]-2))*parms$K[nn]^-parms$n[nn]
               P<-pB[t-agemat,nn]*(PJ/bjoin+cc*(pB[t-agemat,nn]-bjoin)) 
              }
              if(pB[t-agemat,nn]>=bjoin) P<-parms$g[nn]*MSY*(pB[t-agemat,nn]/parms$K[nn])-parms$g[nn]*MSY*(pB[t-agemat,nn]/parms$K[nn])^parms$n[nn] 
           }
          P<-ifelse(P<0,0,P)
          pB[t+1,nn]<-max(0,pB[t,nn]+P-catch1[t-(t-1)])
       } 
    }
  if(grout>=1){     
    plot(apply(pB[star:c(length(pB[,1])),],1,median)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),type="l",
    ylim=c(0,max(apply(pB[star:c(length(pB[,1])),],1,quantile,probs=0.975))*1.20),ylab=paste("Biomass (",grdef$unit,")",sep=""),xlab="Year",
      cex.lab=grdef$cex.lab,
            cex.axis=grdef$cex.axis,main=grdef$mains,cex.main=grdef$cex.main,lwd=grdef$lwd)
    lines(apply(pB[star:c(length(pB[,1])),],1,mean)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=2,lwd=grdef$lwd)
    lines(apply(pB[star:c(length(pB[,1])),],1,quantile,probs=0.025)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=3,lwd=grdef$lwd)
    lines(apply(pB[star:c(length(pB[,1])),],1,quantile,probs=0.975)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=3,lwd=grdef$lwd)
  if(grout==2){
      word.tif = function(filename="Word_Figure_%03d.tif", zoom=tifdef$zoom, 
        width=tifdef$width, height=tifdef$height, pointsize=tifdef$pointsize, ...) {
        if (!grepl("[.]ti[f]+$", filename, ignore.case=TRUE))
        filename = paste0(filename,".tif")
        tiff(filename=filename, compression="lzw", res=96*zoom,
        width=width, height=height, units='cm', pointsize=pointsize, ...)}
     word.tif("ProjBiomass") 
    plot(apply(pB[star:c(length(pB[,1])),],1,median)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),type="l",
    ylim=c(0,max(apply(pB[star:c(length(pB[,1])),],1,quantile,probs=0.975))*1.20),ylab=paste("Biomass (",grdef$unit,")",sep=""),xlab="Year",
      cex.lab=grdef$cex.lab,
            cex.axis=grdef$cex.axis,main=grdef$mains,cex.main=grdef$cex.main,lwd=grdef$lwd)
   lines(apply(pB[star:c(length(pB[,1])),],1,mean)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=2,lwd=grdef$lwd)
     lines(apply(pB[star:c(length(pB[,1])),],1,quantile,probs=0.025)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=3,lwd=grdef$lwd)
     lines(apply(pB[star:c(length(pB[,1])),],1,quantile,probs=0.975)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=3,lwd=grdef$lwd)
   dev.off()
  }
 }
dd<-pB[star:c(length(pB[,1])),]
 rownames(dd)<-c(dlobj$end1yr:c(dlobj$end1yr+projyears))
 ans<-list(dlobj$type,dd)
names(ans)<-c("type","ProjBio")
return(ans)
}# dbsra
#####################################################  
################## Catch msy #######################
if(dlobj$type=="catchmsy"){
     #projtype(0=median MSY, 1=mean MSY, 2=User input) 
      if(projtype<2){ 
       if(projtype==0) catch1<-dlobj$Estimates[1,2] #use median as catch 
       if(projtype==1) catch1<-dlobj$Estimates[1,1] #use mean
       catch1<-rep(catch1,projyears)
      } 
   if(projtype==2){# user input
      if(is.null(projcatch)) stop("projcatch must be specified")
     if(length(projcatch)>1 && length(projcatch)<projyears) stop("length of projcatch does not match length of projyears")
     if(length(projcatch)>1)  catch1<-projcatch
     if(length(projcatch)==1) catch1<-rep(projcatch,projyears)
   }
  
    parms<-dlobj$Values[dlobj$Values[,1]==1,c(3,4)]
    names(parms)<-c("K","r")
    bigB<-read.csv("Biotraj-cmsy.csv",header=FALSE)[,-1]
    bigB<-t(bigB)    
    pB<-as.data.frame(bigB[,c(which(dlobj$Values[,1]==1))])
    pB<-pB[c((length(pB[,1])):length(pB[,1])),]
    pB[c(nrow(pB)+1):c(c(nrow(pB)+1)+projyears-1),]<-NA
   for(nn in 1:nrow(parms)){
     for(t in 1:c(nrow(pB)-1)){ 
     P<-parms$r[nn]*pB[t,nn]*(1-(pB[t,nn]/parms$K[nn]))
     P<-ifelse(P<0,0,P)
     pB[t+1,nn]<-max(0,pB[t,nn]+P-catch1[t])    
    }
  }

if(grout>=1){     
    plot(apply(pB,1,median)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),type="l",
    ylim=c(0,max(apply(pB,1,quantile,probs=0.975))*1.20),ylab=paste("Biomass (",grdef$unit,")",sep=""),xlab="Year",
      cex.lab=grdef$cex.lab,
            cex.axis=grdef$cex.axis,main=grdef$mains,cex.main=grdef$cex.main,lwd=grdef$lwd)
    lines(apply(pB,1,mean)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=2,lwd=grdef$lwd)
    lines(apply(pB,1,quantile,probs=0.025)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=3,lwd=grdef$lwd)
    lines(apply(pB,1,quantile,probs=0.975)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=3,lwd=grdef$lwd)
  if(grout==2){
      word.tif = function(filename="Word_Figure_%03d.tif", zoom=tifdef$zoom, 
        width=tifdef$width, height=tifdef$height, pointsize=tifdef$pointsize, ...) {
        if (!grepl("[.]ti[f]+$", filename, ignore.case=TRUE))
        filename = paste0(filename,".tif")
        tiff(filename=filename, compression="lzw", res=96*zoom,
        width=width, height=height, units='cm', pointsize=pointsize, ...)}
     word.tif("ProjBiomass") 
    plot(apply(pB,1,median)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),type="l",
    ylim=c(0,max(apply(pB,1,quantile,probs=0.975))*1.20),ylab=paste("Biomass (",grdef$unit,")",sep=""),xlab="Year",
      cex.lab=grdef$cex.lab,
            cex.axis=grdef$cex.axis,main=grdef$mains,cex.main=grdef$cex.main,lwd=grdef$lwd)
   lines(apply(pB,1,mean)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=2,lwd=grdef$lwd)
     lines(apply(pB,1,quantile,probs=0.025)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=3,lwd=grdef$lwd)
     lines(apply(pB,1,quantile,probs=0.975)~c(dlobj$end1yr:c(dlobj$end1yr+projyears)),lty=3,lwd=grdef$lwd)
   dev.off()
  }
 }
rownames(pB)<-c(dlobj$end1yr:c(dlobj$end1yr+projyears))
ans<-list(dlobj$type,pB)
names(ans)<-c("type","ProjBio")
 return(ans)  
}
}#function
  