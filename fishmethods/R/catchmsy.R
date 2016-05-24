catchmsy<-function(year = NULL, catch = NULL, catchCV = NULL, catargs = list(dist="none",low=0,up=Inf,unit="MT"),
  l0 = list(low=0,up=1,step=0), lt = list(low=0,up=1,refyr=NULL), sigv = 0,k = list(dist="unif",low=0,up=1,mean=0,sd=0),
  r = list(dist="unif",low=0,up=1,mean=0,sd=0), 
  M = list(dist="unif",low=0.2,up=0.2,mean=0.00,sd=0.00),nsims = 10000, 
  catchout=0,grout = 1,graphs = c(1,2,3,4,5,6,7,8,9,10,11), 
  grargs = list(lwd=1,pch=16,cex=1,nclasses=20,mains=" ",cex.main=1,cex.axis=1,cex.lab=1),
  pstats=list(ol=1,mlty=1,mlwd=1.5,llty=3,llwd=1,ulty=3,ulwd=1),
  grtif=list(zoom=4,width=11,height=13,pointsize=10))
  {
  if(is.null(catch)) stop("No catch data")
  if(length(year)!=length(catch)) stop("Length of year and catch differ")
  #Default setting
  catdef=list(dist="none",low=0,up=Inf,unit="MT")
  if(any(names(catargs)=="dist")){ 
            if(catdef$dist!=catargs$dist) catdef$dist<-catargs$dist
   }
    if(any(names(catargs)=="low")){ 
            if(any(catdef$low!=catargs$low)) catdef$low<-catargs$low
   }
   if(any(names(catargs)=="up")){ 
            if(any(catdef$up!=catargs$up)) catdef$up<-catargs$up
   }
   if(any(names(catargs)=="unit")){ 
            if(catdef$unit!=catargs$unit) catdef$unit<-catargs$unit
   }
  #l0
   l0def=list(low=0,up=1,step=0)
   if(any(names(l0)=="low")){ 
            if(l0def$low!=l0$low) l0def$low<-l0$low
   }
   if(any(names(l0)=="up")){ 
            if(l0def$up!=l0$up) l0def$up<-l0$up
   }
   if(any(names(l0)=="step")){ 
            if(l0def$step!=l0$step) l0def$step<-l0$step
   }
  #lt
  ltdef=list(low=0,up=1,refyr=1)
   if(any(names(lt)=="low")){ 
            if(ltdef$low!=lt$low) ltdef$low<-lt$low
   }
   if(any(names(lt)=="up")){ 
            if(ltdef$up!=lt$up) ltdef$up<-lt$up
   }
   if(any(names(lt)=="refyr")){ 
            if(ltdef$refyr!=lt$refyr) ltdef$refyr<-lt$refyr
   }
  #k
  kdef=list(dist="unif",low=0,up=1,mean=0,sd=0)
   if(any(names(k)=="dist")){ 
            if(kdef$dist!=k$dist) kdef$dist<-k$dist
   }
  if(any(names(k)=="low")){ 
            if(kdef$low!=k$low) kdef$low<-k$low
   }
   if(any(names(k)=="up")){ 
            if(kdef$up!=k$up) kdef$up<-k$up
   }
   if(any(names(k)=="mean")){ 
            if(kdef$mean!=k$mean) kdef$mean<-k$mean
   }
   if(any(names(k)=="sd")){ 
            if(kdef$sd!=k$sd) kdef$sd<-k$sd
   }
  #r
  rdef=list(dist="unif",low=0,up=1,mean=0,sd=0)
  if(any(names(r)=="dist")){ 
            if(rdef$dist!=r$dist) rdef$dist<-r$dist
   }
  if(any(names(r)=="low")){ 
            if(rdef$low!=r$low) rdef$low<-r$low
   }
   if(any(names(r)=="up")){ 
            if(rdef$up!=r$up) rdef$up<-r$up
   }
   if(any(names(r)=="mean")){ 
            if(rdef$mean!=r$mean) rdef$mean<-r$mean
   }
   if(any(names(r)=="sd")){ 
            if(rdef$sd!=r$sd) rdef$sd<-r$sd
   }
 
 #M
  Mdef=list(dist="unif",low=0.2,up=0.2,mean=0.00,sd=0.00)
   if(any(names(M)=="dist")){ 
            if(Mdef$dist!=M$dist) Mdef$dist<-M$dist
   }
  if(any(names(M)=="low")){ 
            if(Mdef$low!=M$low) Mdef$low<-M$low
   }
   if(any(names(M)=="up")){ 
            if(Mdef$up!=M$up) Mdef$up<-M$up
   }
   if(any(names(M)=="mean")){ 
            if(Mdef$mean!=M$mean) Mdef$mean<-M$mean
   }
   if(any(names(M)=="sd")){ 
            if(Mdef$sd!=M$sd) Mdef$sd<-M$sd
   }
  #gr def
  grdef=list(pch=16,cex.axis=1,cex.lab=1,cex=1,nclasses=20,
        mains=" ",cex.main=1,lwd=1)
   if(any(names(grargs)=="pch")){ 
            if(any(grdef$pch!=grargs$pch)) grdef$pch<-grargs$pch
   }
  if(any(names(grargs)=="cex.axis")){ 
            if(any(grdef$cex.axis!=grargs$cex.axis)) grdef$cex.axis<-grargs$cex.axis
   }
   if(any(names(grargs)=="cex.lab")){ 
            if(any(grdef$cex.lab!=grargs$cex.lab)) grdef$cex.lab<-grargs$cex.lab
   }
   if(any(names(grargs)=="cex")){ 
            if(any(grdef$cex!=grargs$cex)) grdef$cex<-grargs$cex
   }
   if(any(names(grargs)=="nclasses")){ 
            if(any(grdef$nclasses!=grargs$nclasses)) grdef$nclasses<-grargs$nclasses
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
  #pstat
  pstdef=list(ol=1,mlty=1,mlwd=1.5,llty=3,llwd=1,ulty=3,ulwd=1)
   if(any(names(pstats)=="ol")){ 
            if(pstdef$ol!=pstats$ol) pstdef$ol<-pstats$ol
   }
  if(any(names(pstats)=="mlty")){ 
            if(pstdef$mlty!=pstats$mlty) pstdef$mlty<-pstats$mlty
   }
   if(any(names(pstats)=="mlwd")){ 
            if(pstdef$mlwd!=pstats$mlwd) pstdef$mlwd<-pstats$mlwd
   }
   if(any(names(pstats)=="llty")){ 
            if(pstdef$llty!=pstats$llty) pstdef$llty<-pstats$llty
   }
   if(any(names(pstats)=="llwd")){ 
            if(pstdef$llwd!=pstats$llwd) pstdef$llwd<-pstats$llwd
   }
   if(any(names(pstats)=="ulty")){ 
            if(pstdef$ulty!=pstats$ulty) pstdef$ulty<-pstats$ulty
   }
   if(any(names(pstats)=="ulwd")){ 
            if(pstdef$ulwd!=pstats$ulwd) pstdef$ulwd<-pstats$ulwd
   }
  #tif arguments
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

  if(l0def$low<0|l0def$low>1) stop("l0 low can only range from 0 to 1")
  if(l0def$up<0|l0def$up>1) stop("l0 up can only range from 0 to 1")
  if(l0def$step<0) stop("l0 step must be >=0.")
  if(ltdef$low<0|ltdef$low>1) stop("lt low can only range from 0 to 1")
  if(ltdef$up<0|ltdef$up>1) stop("lt up can only range from 0 to 1")
 if(any(is.na(catch))) stop("Missing catch values are not allowed.")
  
  if(catdef$dist!="none"){
    if(catdef$dist %in% c("norm","lnorm") & any(!is.numeric(catchCV))) stop("catchCV is required for resampling")
    catdata<-as.data.frame(cbind(year,catch,catchCV))
  }
  if(catdef$dist=="none"){
   catdata<-as.data.frame(cbind(year,catch))
  }    
  word.tif = function(filename="Word_Figure_%03d.tif", zoom=tifdef$zoom, 
    width=tifdef$width, height=tifdef$height, pointsize=tifdef$pointsize, ...) {
    if (!grepl("[.]ti[f]+$", filename, ignore.case=TRUE))
      filename = paste0(filename,".tif")
       tiff(filename=filename, compression="lzw", res=96*zoom,
       width=width, height=height, units='cm', pointsize=pointsize, ...)}
#random ditributions
  getrandom<-function(n,spec,a=-Inf,b=Inf,mean=NULL,sd=NULL){
    if(!spec %in% c("unif","beta","lnorm","norm","gamma")) stop ("Unknown distribution name.")
    if(spec %in% c("beta","lnorm","norm","gamma")){
      if(is.null(mean)|is.null(sd)) stop("mean and sd must be specified.")
       if(is.na(mean)|is.na(sd)) stop("mean and sd must be specified.")
    }
  if(spec=="beta"){
    if(is.null(a)|is.null(b)|a==-Inf|b==Inf) stop("lower and upper limits are required for the beta distribution.")
    if(is.na(a)|is.na(b)) stop("lower and upper limits are required for the beta distribution.")
    }  
  if(spec=="unif"){
    if(is.na(a)|is.na(b)|is.null(a)|is.null(b)) stop("unif requires specified lower and upper values")
    if(a==-Inf|b==Inf) stop("unif requires specified lower and upper values")
  }
  
qtrunc<-function(p,spec,a=-Inf,b=Inf,...){
    ttt<-p
    G<-get(paste("p",spec,sep=""),mode="function")
    Gin<-get(paste("q",spec,sep=""),mode="function")
    ttt<-Gin(G(a,...)+p*(G(b,...)-G(a,...)),...)
    return(ttt)
  }

  rtrunc<-function(n,spec,a=-Inf,b=Inf,...){
    x<-u<-runif(n,min=0,max=1)
    x<-qtrunc(u,spec,a=a,b=b,...)
    return(x)
  }
if(spec=="beta"){
    return(rtrunc(n,"beta",a=a,b=b,shape1=mean*(((mean*(1-mean))/sd^2)-1),
      shape2=(1-mean)*(((mean*(1-mean))/sd^2)-1)))
  }
  if(spec=="unif") return(rtrunc(n,"unif",min=a,max=b))
  if(spec=="norm") return(rtrunc(n,"norm",a=a,b=b,mean=mean,sd=sd))
  if(spec=="lnorm") return(rtrunc(n,"lnorm",a=a,b=b,meanlog=mean,sdlog=sd))
  if(spec=="gamma") return(rtrunc(n,"gamma",a=a,b=b,shape=(mean/sd)^2,scale=sd^2/mean))
} 
   
   timelen<-length(year)
   termyear<-NULL
    if(l0def$step>0){
      steps<-seq(l0def$low,l0def$up,l0def$step) 
      storep<-matrix(0,nrow=nsims*c(length(steps)),ncol=11) #l0, k,r,BK,M,likel,MSY,Fmsy,Umsy,OFL,Bt/k
    }
    if(l0def$step==0) storep<-matrix(0,nrow=nsims,ncol=11) #l0, k,r,BK,M,likel,MSY,Fmsy,Umsy,OFL,Bt/k
     termyear<-length(year)+1
     refyr<-which(lt$refyr==c(year,year[length(year)]+1))
     kk<-kdef$mean
     rr<-rdef$mean
     MM<-Mdef$mean
     lt$refyr<-1999
  #Program
    cnt<-0 
    for(n in 1:nsims){
     B<-NULL
      if(kdef$dist!="none") kk<-getrandom(1,kdef$dist,a=kdef$low,b=kdef$up,mean=kdef$mean,sd=kdef$sd)
      if(rdef$dist!="none") rr<-getrandom(1,rdef$dist,a=rdef$low,b=rdef$up,mean=rdef$mean,sd=rdef$sd)
      if(Mdef$dist!="none") MM<-getrandom(1,Mdef$dist,a=Mdef$low,b=Mdef$up,mean=Mdef$mean,sd=Mdef$sd)
      ## Resample catch if desired    
       dcatch<-NULL
      if(catdef$dist=="none"){
        dcatch<-catdata[,2] 
        catchout<-0
      }      
      if(catdef$dist!="none"){
        for(cc in 1:c(length(catdata[,2]))){
          if(catdef$dist=="unif") dcatch[cc]<-getrandom(1,catdef$dist,a=catdef$low[cc],b=catdef$up[cc])
          if(catdef$dist=="norm") dcatch[cc]<-getrandom(1,catdef$dist,a=catdef$low,b=catdef$up,mean=catdata[cc,2],sd=catdata[cc,2]*catdata[cc,3])
          if(catdef$dist=="lnorm") dcatch[cc]<-getrandom(1,catdef$dist,a=log(catdef$low),b=log(catdef$up),mean=log(catdata[cc,2]),sd=sqrt(log(catdata[cc,3]^2+1)))
        }
        dcatch<-ifelse(dcatch<0|is.infinite(dcatch),0,dcatch)
       
      }
          
     if(l0def$step==0) steps<-runif(1,min=l0def$low,max=l0def$up)     
  
     for(sbio in 1:c(length(steps))){# loop throgh l0 values
       cnt<-cnt+1
       B[1]<-steps[sbio]*kk*exp(rnorm(1,0,sigv))
       for(t in 1:timelen) B[t+1]<-max(0,(B[t]+rr*B[t]*(1-(B[t]/kk))-dcatch[t])*exp(rnorm(1,0,sigv))) 
       bll<-0
       if((B[refyr]/kk)>=ltdef$low && (B[refyr]/kk)<=ltdef$up && min(B)>0 && max(B)<=kk) bll<-1 
       if(n==1){
         write.table(t(c(bll,B)),file="Biotraj-cmsy.csv",sep=",",row.names=FALSE,col.names=FALSE,append=FALSE)
         if(catchout==1) write.table(t(c(bll,dcatch)),file="Catchtraj-cmsy.csv",sep=",",row.names=FALSE,col.names=FALSE,append=FALSE)
       }  
         
       if(n>1){
         write.table(t(c(bll,B)),file="Biotraj-cmsy.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
         if(catchout==1)  write.table(t(c(bll,dcatch)),file="Catchtraj-cmsy.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
       }      
       storep[cnt,1]<-bll
        storep[cnt,2]<-steps[sbio]
        storep[cnt,3]<-kk
        storep[cnt,4]<-rr
        storep[cnt,5]<-MM
        storep[cnt,6]<-(rr*kk/4) #MSY
        storep[cnt,7]<-kk/2 #Bmsy
        storep[cnt,8]<-rr/2 #Fmsy
        storep[cnt,9]<-storep[cnt,8]/(storep[cnt,8]+MM)*(1-exp(-MM-storep[cnt,8]))#Umsy
        storep[cnt,10]<-B[termyear]*storep[cnt,9] #OFL
        storep[cnt,11]<-B[refyr]/kk
    }#steps
  }#n
   # Get results for bll==1 for graphing
  storep<-as.data.frame(storep)
  colnames(storep)<-c("likeli","l0","k","r","M","MSY","Bmsy","Fmsy","Umsy","OFL","Btk")
  datar<-storep[storep[,1]==1,]
   
  if(length(datar[,1])>0){ 
    mMSY<-round(mean(datar$MSY),3)
    MSY95<-round(quantile(datar$MSY,probs=c(0.025,0.50,0.975)),4)
    mk<-round(mean(datar$k),3)
    k95<-round(quantile(datar$k,probs=c(0.025,0.50,0.975)),4)
    mr<-round(mean(datar$r),3)
    r95<-round(quantile(datar$r,probs=c(0.025,0.50,0.975)),4)
    mBMSY<-round(mean(datar$Bmsy),3)
    BMSY95<-round(quantile(datar$Bmsy,probs=c(0.025,0.50,0.975)),4)
    mFMSY<-round(mean(datar$Fmsy),3)
    FMSY95<-round(quantile(datar$Fmsy,probs=c(0.025,0.50,0.975)),4)
    mUMSY<-round(mean(datar$Umsy),3)
    UMSY95<-round(quantile(datar$Umsy,probs=c(0.025,0.50,0.975)),4)
    mOFL<-round(mean(datar$OFL),3)
    OFL95<-round(quantile(datar$OFL,probs=c(0.025,0.50,0.975)),4)
     mM<-round(mean(datar$M),3)
    M95<-round(quantile(datar$M,probs=c(0.025,0.50,0.975)),4)
    
    ### Graphs######################
   
    # Plot graphs
  if(grout>0){   
     grunits<-data.frame(gr=c(1:11),cexa=0,cexl=0,cexx=0,nclass=0,mains=" ",cexmain=0,lwd=0,pch=0,stringsAsFactors=FALSE)
     grunits$lwd<-ifelse(grunits$gr %in% c(1,11),grdef$lwd,0) 
     grunits$pch<-ifelse(grunits$gr %in% c(2),grdef$pch,0)   
     grunits[,2]<-grdef$cex.axis
     grunits[,3]<-grdef$cex.lab
     grunits[,4]<-grdef$cex
     grunits$nclass<-ifelse(grunits$gr %in% c(3,4,5,6,7,8,9,10),grdef$nclasses,grunits$nclass)
     grunits$cexmain<-grdef$cex.main
     grunits[graphs,6]<-grdef$mains   
    if(any(graphs==1)){
          plot(catch~year,type="l",xlab="Year",ylab=paste("Catch (",catdef$unit,")",sep=""),
          ylim=c(0,round(max(catch,mMSY,MSY95[3]),0)),cex=grunits[1,4],cex.lab=grunits[1,3],
            cex.axis=grunits[1,2],main=grunits[1,6],cex.main=grunits[1,7],lwd=grunits[1,9])
           if(pstdef$ol==1){
             abline(h=mMSY,lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(h=MSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
             abline(h=MSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
           }
          if(grout==2){
           word.tif("catch-cmsy")
            plot(catch~year,type="l",xlab="Year",ylab=paste("Catch (",catdef$unit,")",sep=""),
            ylim=c(0,round(max(catch,mMSY,MSY95[3]),0)),cex=grunits[1,4],cex.lab=grunits[1,3],
              cex.axis=grunits[1,2],main=grunits[1,6],cex.main=grunits[1,7],lwd=grunits[1,9])
             if(pstdef$ol==1){
             abline(h=mMSY,lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(h=MSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
             abline(h=MSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
             }
           dev.off()
         }
      }
      if(any(graphs==2)){
           plot(datar$k~datar$r,xlim=c(0,max(datar$r)*1.50),
             ylim=c(0,max(datar$k)),ylab="k",xlab="r",col="black",cex=grunits[2,4],cex.lab=grunits[2,3],cex.axis=grunits[2,2],main=grunits[2,6],cex.main=grunits[2,7],pch=grunits[2,8])
          if(grout==2){
            word.tif("kvsr-cmsy")  
            plot(datar$k~datar$r,xlim=c(0,max(datar$r)*1.50),
            ylim=c(0,max(datar$k)),ylab="k",xlab="r",col="black",cex=grunits[2,4],cex.lab=grunits[2,3],cex.axis=grunits[2,2],main=grunits[2,6],cex.main=grunits[2,7],pch=grunits[2,8])
            dev.off()
          }
      }  
      if(any(graphs==3)){
        hist(datar$r,freq=FALSE,xlim=c(0,max(datar$r)*1.50),xlab="r",nclass=grunits[3,5],cex.lab=grunits[3,3],cex.axis=grunits[3,2],main=grunits[3,6],cex.main=grunits[3,7])
           if(pstdef$ol==1){ 
             abline(v=mr,lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(v=r95[1],lwd=pstdef$lwd,lty=pstdef$llty)
             abline(v=r95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
            }
        if(grout==2){
            word.tif("rden-cmsy") 
             hist(datar$r,freq=FALSE,xlim=c(0,max(datar$r)*1.50),xlab="r",nclass=grunits[3,5],cex.lab=grunits[3,3],cex.axis=grunits[3,2],main=grunits[3,6],cex.main=grunits[3,7])
            if(pstdef$ol==1){ 
             abline(v=mr,lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(v=r95[1],lwd=pstdef$llwd,lty=pstdef$llty)
             abline(v=r95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
            }
          dev.off()
         }
      }
      if(any(graphs==4)){
           hist(datar$k,freq=FALSE,xlim=c(0,max(datar$k)*1.50),xlab="k",nclass=grunits[4,5],cex.lab=grunits[4,3],cex.axis=grunits[4,2],main=grunits[4,6],cex.main=grunits[4,7])
            if(pstdef$ol==1){ 
             abline(v=mk,lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(v=k95[1],lwd=pstdef$llwd,lty=pstdef$llty)
             abline(v=k95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
            }
        if(grout==2){
           word.tif("kden-cmsy") 
           hist(datar$k,freq=FALSE,xlim=c(0,max(datar$k)*1.50),xlab="k",nclass=grunits[4,5],cex.lab=grunits[4,3],cex.axis=grunits[4,2],main=grunits[4,6],cex.main=grunits[4,7])
            if(pstdef$ol==1){ 
             abline(v=mk,lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(v=k95[1],lwd=pstdef$llwd,lty=pstdef$llty)
             abline(v=k95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
            }
          dev.off()
         }
      }    
     
    if(any(graphs==5)){      
          hist(datar$M,freq=FALSE,xlim=c(min(datar$M)*0.80,max(datar$M)*1.20),xlab="M",nclass=grunits[5,5],cex.lab=grunits[5,3],cex.axis=grunits[5,2],main=grunits[5,6],cex.main=grunits[5,7])
          if(pstdef$ol==1){ 
           abline(v=mM,lwd=pstdef$mlwd,lty=pstdef$mlty)
           abline(v=M95[1],lwd=pstdef$llwd,lty=pstdef$llty)
           abline(v=M95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
          }
         if(grout==2){
           word.tif("Mden-cmsy") 
            hist(datar$M,freq=FALSE,xlim=c(min(datar$M)*0.80,max(datar$M)*1.20),xlab="M",nclass=grunits[5,5],cex.lab=grunits[5,3],cex.axis=grunits[5,2],main=grunits[5,6],cex.main=grunits[5,7])
            if(pstdef$ol==1){ 
             abline(v=mM,lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(v=M95[1],lwd=pstdef$llwd,lty=pstdef$llty)
             abline(v=M95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
            }
           dev.off()
          }
      }  
     if(any(graphs==6)){      
         hist(datar$MSY,freq=FALSE,xlim=c(min(datar$MSY)*0.80,max(datar$MSY)*1.20),xlab=paste("MSY (",catargs[4],")",sep=""),nclass=grunits[6,5],cex.lab=grunits[6,3],cex.axis=grunits[6,2],main=grunits[6,6],cex.main=grunits[6,7])
         if(pstdef$ol==1){ 
           abline(v=mMSY,lwd=pstdef$mlwd,lty=pstdef$mlty)
           abline(v=MSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
           abline(v=MSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
        if(grout==2){
           word.tif("MSYden-cmsy") 
             hist(datar$MSY,freq=FALSE,xlim=c(min(datar$MSY)*0.80,max(datar$MSY)*1.20),xlab=paste("MSY (",catargs[4],")",sep=""),nclass=grunits[6,5],cex.lab=grunits[6,3],cex.axis=grunits[6,2],main=grunits[6,6],cex.main=grunits[6,7])
              if(pstdef$ol==1){ 
               abline(v=mMSY,lwd=pstdef$mlwd,lty=pstdef$mlty)
               abline(v=MSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
               abline(v=MSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
              }
           dev.off()
          }
     }
    if(any(graphs==7)){   
       hist(datar$Bmsy,freq=FALSE,xlim=c(0,max(datar$Bmsy)*1.50),xlab=paste("Bmsy (",catdef$unit,")",sep=""),nclass=grunits[7,5],cex.lab=grunits[7,3],cex.axis=grunits[7,2],main=grunits[7,6],cex.main=grunits[7,7])
        if(pstdef$ol==1){
          abline(v=mBMSY,lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=BMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=BMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
        }
      if(grout==2){
           word.tif("Bmsyden-cmsy") 
               hist(datar$Bmsy,freq=FALSE,xlim=c(0,max(datar$Bmsy)*1.50),xlab=paste("Bmsy (",catdef$unit,")",sep=""),nclass=grunits[7,5],cex.lab=grunits[7,3],cex.axis=grunits[7,2],main=grunits[7,6],cex.main=grunits[7,7])
           if(pstdef$ol==1){ 
            abline(v=mBMSY,lwd=pstdef$mlwd,lty=pstdef$mlty)
            abline(v=BMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
            abline(v=BMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
           }
           dev.off()
      }
    }  
    if(any(graphs==8)){   
         hist(datar$Fmsy,freq=FALSE,xlim=c(0,max(datar$Fmsy)*1.50),xlab="Fmsy",nclass=grunits[8,5],cex.lab=grunits[8,3],cex.axis=grunits[8,2],main=grunits[8,6],cex.main=grunits[8,7])
         if(pstdef$ol==1){
          abline(v=mFMSY,lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=FMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=FMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
        if(grout==2){
           word.tif("Fmsyden-cmsy") 
           hist(datar$Fmsy,freq=FALSE,xlim=c(0,max(datar$Fmsy)*1.50),xlab="Fmsy",nclass=grunits[8,5],cex.lab=grunits[8,3],cex.axis=grunits[8,2],main=grunits[8,6],cex.main=grunits[8,7])
            if(pstdef$ol==1){
            abline(v=mFMSY,lwd=pstdef$mlwd,lty=pstdef$mlty)
            abline(v=FMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
            abline(v=FMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
           }
          dev.off()
        }
    }
    if(any(graphs==9)){   
         hist(datar$Umsy,freq=FALSE,xlim=c(0,max(datar$Umsy)*1.50),xlab="Umsy",nclass=grunits[9,5],cex.lab=grunits[9,3],cex.axis=grunits[9,2],main=grunits[9,6],cex.main=grunits[9,7])
         if(pstdef$ol==1){
          abline(v=mUMSY,lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=UMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=UMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
        if(grout==2){
           word.tif("Umsyden-cmsy") 
           hist(datar$Umsy,freq=FALSE,xlim=c(0,max(datar$Umsy)*1.50),xlab="Umsy",nclass=grunits[9,5],cex.lab=grunits[9,3],cex.axis=grunits[9,2],main=grunits[9,6],cex.main=grunits[9,7])
            if(pstdef$ol==1){
            abline(v=mUMSY,lwd=pstdef$mlwd,lty=pstdef$mlty)
            abline(v=UMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
            abline(v=UMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
           }
          dev.off()
        }
    }
       if(any(graphs==10)){   
         hist(datar$OFL,freq=FALSE,xlim=c(0,max(datar$OFL)*1.50),xlab=paste("OFL (",catdef$unit,")",sep=""),nclass=grunits[10,5],cex.lab=grunits[10,3],cex.axis=grunits[10,2],main=grunits[10,6],cex.main=grunits[10,7])
         if(pstdef$ol==1){
          abline(v=mOFL,lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=OFL95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=OFL95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
        if(grout==2){
           word.tif("OFLden-cmsy") 
           hist(datar$OFL,freq=FALSE,xlim=c(0,max(datar$OFL)*1.50),xlab=paste("OFL (",catdef$unit,")",sep=""),nclass=grunits[10,5],cex.lab=grunits[10,3],cex.axis=grunits[10,2],main=grunits[10,6],cex.main=grunits[10,7])
            if(pstdef$ol==1){
            abline(v=mOFL,lwd=pstdef$mlwd,lty=pstdef$mlty)
            abline(v=OFL95[1],lwd=pstdef$llwd,lty=pstdef$llty)
            abline(v=OFL95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
           }
          dev.off()
        }
    }
    if(any(graphs==11)){
         bigB<-read.csv("Biotraj-cmsy.csv",header=FALSE)
         ar<-bigB[,1]
         bigB<-t(bigB[,-1])
         cols<-ifelse(ar==0,"gray85","black")
         types<-ifelse(ar==0,3,1)
                 
          par(mfrow=c(1,2))
         matplot(y=bigB[,c(which(ar==1))],x=c(year,year[length(year)]+1),type="l",lty=types[c(which(ar==1))],xlab="Year",ylab=paste("Biomass (",catdef$unit,")",sep=""),
            ylim=c(0,max(bigB)),cex=grunits[11,4],cex.lab=grunits[11,3],col=cols[c(which(ar==1))],main="Accepted",
              cex.axis=grunits[11,2],cex.main=grunits[11,7],lwd=grunits[11,8])    
           lines(y=apply(bigB[,c(which(ar==1))],1,median),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$mlty,col="red")
           lines(y=apply(bigB[,c(which(ar==1))],1,function(x){quantile(x,probs=0.025)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$llty,col="red")
           lines(y=apply(bigB[,c(which(ar==1))],1,function(x){quantile(x,probs=0.975)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$ulty,col="red")
 
         matplot(y=bigB[,c(which(ar==0))],x=c(year,year[length(year)]+1),type="l",lty=types[c(which(ar==0))],xlab="Year",ylab=paste("Biomass (",catdef$unit,")",sep=""),
            ylim=c(0,max(bigB)),cex=grunits[11,4],cex.lab=grunits[11,3],col=cols[c(which(ar==0))],main="Rejected",
              cex.axis=grunits[11,2],cex.main=grunits[11,7],lwd=grunits[11,8])    
           lines(y=apply(bigB[,c(which(ar==0))],1,median),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$mlty,col="red")
           lines(y=apply(bigB[,c(which(ar==0))],1,function(x){quantile(x,probs=0.025)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$llty,col="red")
           lines(y=apply(bigB[,c(which(ar==0))],1,function(x){quantile(x,probs=0.975)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$ulty,col="red")
        if(grout==2){
            word.tif("Biotraj-cmsy") 
            par(mfrow=c(1,2))
            matplot(y=bigB[,c(which(ar==1))],x=c(year,year[length(year)]+1),type="l",lty=types[c(which(ar==1))],xlab="Year",ylab=paste("Biomass (",catdef$unit,")",sep=""),
            ylim=c(0,max(bigB)),cex=grunits[11,4],cex.lab=grunits[11,3],col=cols[c(which(ar==1))],main="Accepted",
              cex.axis=grunits[11,2],cex.main=grunits[11,7],lwd=grunits[11,8])    
           lines(y=apply(bigB[,c(which(ar==1))],1,median),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$mlty,col="red")
           lines(y=apply(bigB[,c(which(ar==1))],1,function(x){quantile(x,probs=0.025)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$llty,col="red")
           lines(y=apply(bigB[,c(which(ar==1))],1,function(x){quantile(x,probs=0.975)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$ulty,col="red")
 
         matplot(y=bigB[,c(which(ar==0))],x=c(year,year[length(year)]+1),type="l",lty=types[c(which(ar==0))],xlab="Year",ylab=paste("Biomass (",catdef$unit,")",sep=""),
            ylim=c(0,max(bigB)),cex=grunits[11,4],cex.lab=grunits[11,3],col=cols[c(which(ar==0))],main="Rejected",
              cex.axis=grunits[11,2],cex.main=grunits[11,7],lwd=grunits[11,8])    
           lines(y=apply(bigB[,c(which(ar==0))],1,median),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$mlty,col="red")
           lines(y=apply(bigB[,c(which(ar==0))],1,function(x){quantile(x,probs=0.025)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$llty,col="red")
           lines(y=apply(bigB[,c(which(ar==0))],1,function(x){quantile(x,probs=0.975)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$ulty,col="red")
          dev.off()
        }
       rm(bigB)
       par(mfrow=c(1,1))
    }
  }#grout>0
    
    outs<-data.frame(Distr=NA,Min=NA,Max=NA,Mean=NA,sd=NA)
    if(l0def$step==0) outs[1,]<-cbind("increment",l0def$low,l0def$up,NA,NA)
    if(l0def$step>0) outs[1,]<-cbind("unif",l0def$low,l0def$up,NA,NA)
    outs[2,]<-cbind(kdef$dist,kdef$low,kdef$up,kdef$mean,kdef$sd)
    outs[3,]<-cbind(rdef$dist,rdef$low,rdef$up,rdef$mean,rdef$sd)
    outs[4,]<-cbind(Mdef$dist,Mdef$low,Mdef$up,Mdef$mean,Mdef$sd)
    outs[5,]<-cbind(NA,ltdef$low,ltdef$up,NA,NA)
    outs[6,]<-cbind(NA,ltdef$refyr,NA,NA,NA)
    colnames(outs)<-c("Distr","Lower","Upper","Mean","SD")
    rownames(outs)<-c("l0","k","r","M","lt","refyr")
  
    outs1<-data.frame(Mean=NA,Median=NA,per2_5=NA,per97_5=NA)
    outs1[1,]<-cbind(mk,k95[[2]],k95[[1]],k95[[3]])
    outs1[2,]<-cbind(mr,r95[[2]],r95[[1]],r95[[3]])
    outs1[3,]<-cbind(mM,M95[[2]],M95[[1]],M95[[3]])
    colnames(outs1)<-c("Mean (ll=1)","Median (ll=1)","2.5% (ll=1)","97.5% (ll=1)")
    rownames(outs1)<-c("k","r","M")
    
    outs2<-data.frame(Mean=NA,Median=NA, per2_5=NA,per97_5=NA)
    outs2[1,]<-cbind(mMSY,MSY95[[2]],MSY95[[1]],MSY95[[3]])
    outs2[2,]<-cbind(mBMSY,BMSY95[[2]],BMSY95[[1]],BMSY95[[3]])
    outs2[3,]<-cbind(mFMSY,FMSY95[[2]],FMSY95[[1]],FMSY95[[3]])
    outs2[4,]<-cbind(mUMSY,UMSY95[[2]],UMSY95[[1]],UMSY95[[3]])
    outs2[5,]<-cbind(mOFL,OFL95[[2]],OFL95[[1]],OFL95[[3]])
    colnames(outs2)<-c("Mean","Median","2.5%","97.5%")
    rownames(outs2)<-c("MSY","BMSY","FMSY","Umsy","OFL")
    ans<-list(outs,outs1,outs2,storep,max(year)+1,"catchmsy")
    names(ans)<-c("Initial","Parameters","Estimates","Values","end1yr","type")
    
} #length(datar[,1])>1
  
if(length(datar[,1])==0){
      outs<-data.frame(Distr=NA,Min=NA,Max=NA,Mean=NA,sd=NA)
    if(l0def$step==0) outs[1,]<-cbind("increment",l0def$low,l0def$up,NA,NA)
    if(l0def$step>0) outs[1,]<-cbind("unif",l0def$low,l0def$up,NA,NA)
    outs[2,]<-cbind(kdef$dist,kdef$low,kdef$up,kdef$mean,kdef$sd)
    outs[3,]<-cbind(rdef$dist,rdef$low,rdef$up,rdef$mean,rdef$sd)
    outs[4,]<-cbind(Mdef$dist,Mdef$low,Mdef$up,Mdef$mean,Mdef$sd)
    outs[5,]<-cbind(NA,ltdef$low,ltdef$up,NA,NA)
    outs[6,]<-cbind(NA,ltdef$refyr,NA,NA,NA)
    colnames(outs)<-c("Distr","Lower","Upper","Mean","SD")
    rownames(outs)<-c("l0","k","r","M","lt","refyr")
    outs1<-data.frame(Mean=NA,Median=NA,per2_5=NA,per97_5=NA)
    colnames(outs1)<-c("Mean (ll=1)","Median (ll=1)","2.5% (ll=1)","97.5% (ll=1)")
    outs2<-data.frame(Mean=NA,Median=NA,per2_5=NA,per97_5=NA)
    colnames(outs2)<-c("Mean","Median","2.5%","97.5%")
    ans<-list(outs,outs1,outs2,storep,max(year)+1,"catchmsy")
    names(ans)<-c("Initial","Parameters","Estimates","Values","end1yr","type")
     warning("None of the runs had a likelihood equal to 1")
  }
    return(ans) 
}#function
