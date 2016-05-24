dbsra<-function(year = NULL, catch = NULL, catchCV = NULL, 
  catargs = list(dist="none",low=0,up=Inf,unit="MT"),
  agemat=NULL,
  k = list(low=0,up=NULL,tol=0.01,permax=1000), 
  btk = list(dist="unif",low=0,up=1,mean=0,sd=0,refyr=NULL),
  fmsym = list(dist="unif",low=0,up=1,mean=0,sd=0),
  bmsyk = list(dist="unif",low=0,up=1,mean=0,sd=0),
  M = list(dist="unif",low=0.0,up=1,mean=0,sd=0),
  nsims = 10000, catchout=0,
  grout = 1,
  graphs = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14), 
  grargs = list(lwd=1,cex=1,nclasses=20,mains=" ",cex.main=1,cex.axis=1,cex.lab=1),
  pstats=list(ol=1,mlty=1,mlwd=1.5,llty=3,llwd=1,ulty=3,ulwd=1),
  grtif=list(zoom=4,width=11,height=13,pointsize=10))
  {
  if(is.null(catch)) stop("No catch data")
  if(length(year)!=length(catch)) stop("Length of year and catch differ")
    if(btk$refyr>c(max(year)+1)) stop("refyr is beyond the max year allowed (i.e.,max(year)+1)") 
  #Default setting
  catdef=list(dist="none",low=0,up=Inf,unit="MT")
  if(any(is.na(catch))) stop("There are missing values in catch")
  if(any(is.na(year))) stop("There are missing values in year")
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
  #kintial
   kdef=list(low=0,up=max(catch),tol=1,permax=1000)
   if(any(names(k)=="low")){ 
            if(kdef$low!=k$low) kdef$low<-k$low
   }
   if(any(names(k)=="up")){ 
            if(kdef$up!=k$up) kdef$up<-k$up
   }
   if(any(names(k)=="tol")){ 
            if(kdef$tol!=k$tol) kdef$tol<-k$tol
   }
 if(any(names(k)=="permax")){ 
            if(kdef$permax!=k$permax) kdef$permax<-k$permax
   }
  #Fmsy over M
  fmsymdef=list(dist="unif",low=0,up=1,mean=0,sd=0)
   if(any(names(fmsym)=="dist")){ 
            if(fmsymdef$dist!=fmsym$dist) fmsymdef$dist<-fmsym$dist
   }
  if(any(names(fmsym)=="low")){ 
            if(fmsymdef$low!=fmsym$low) fmsymdef$low<-fmsym$low
   }
   if(any(names(fmsym)=="up")){ 
            if(fmsymdef$up!=fmsym$up) fmsymdef$up<-fmsym$up
   }
   if(any(names(fmsym)=="mean")){ 
            if(fmsymdef$mean!=fmsym$mean) fmsymdef$mean<-fmsym$mean
   }
   if(any(names(fmsym)=="sd")){ 
            if(fmsymdef$sd!=fmsym$sd) fmsymdef$sd<-fmsym$sd
   }
  #btk
  btkdef=list(dist="unif",low=0,up=1,mean=0,sd=0,refyr=max(year))
  if(any(names(btk)=="dist")){ 
            if(btkdef$dist!=btk$dist) btkdef$dist<-btk$dist
   }
  if(any(names(btk)=="low")){ 
            if(btkdef$low!=btk$low) btkdef$low<-btk$low
   }
   if(any(names(btk)=="up")){ 
            if(btkdef$up!=btk$up) btkdef$up<-btk$up
   }
   if(any(names(btk)=="mean")){ 
            if(btkdef$mean!=btk$mean) btkdef$mean<-btk$mean
   }
   if(any(names(btk)=="sd")){ 
            if(btkdef$sd!=btk$sd) btkdef$sd<-btk$sd
   }
   if(any(names(btk)=="refyr")){ 
            if(btkdef$refyr!=btk$refyr) btkdef$refyr<-btk$refyr
   }
  
  bmsykdef=list(dist="unif",low=0.5,up=0.5,mean=0,sd=0)
   if(any(names(bmsyk)=="dist")){ 
            if(bmsykdef$dist!=bmsyk$dist) bmsykdef$dist<-bmsyk$dist
   }
  if(any(names(bmsyk)=="low")){ 
            if(bmsykdef$low!=bmsyk$low) bmsykdef$low<-bmsyk$low
   }
   if(any(names(bmsyk)=="up")){ 
            if(bmsykdef$up!=bmsyk$up) bmsykdef$up<-bmsyk$up
   }
   if(any(names(bmsyk)=="mean")){ 
            if(bmsykdef$mean!=bmsyk$mean) bmsykdef$mean<-bmsyk$mean
   }
   if(any(names(bmsyk)=="sd")){ 
            if(bmsykdef$sd!=bmsyk$sd) bmsykdef$sd<-bmsyk$sd
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
  grdef=list(lwd=1,cex.axis=1,cex.lab=1,cex=1,nclasses=20,
        mains=" ",cex.main=1)
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

  if(btkdef$low<0|btkdef$low>1) stop("btk low can only range from 0 to 1")
  if(btkdef$up<0|btkdef$up>1) stop("btk up can only range from 0 to 1")
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
   if(catdef$dist=="unif"){
         if(length(catdef$low)!=length(catdata[,2])|
            length(catdef$up)!=length(catdata[,2])) stop("The length of catargs$low and/or catargs$up should be the same length as catch")
    }   
   timelen<-length(year)
   refyr<-which(btk$refyr==c(year,year[length(year)]+1))
   
  storep<-data.frame(NULL)
  ###########Program
  fmsymX<-fmsymdef$mean
  btkX<-btkdef$mean 
  bmsykX<-bmsykdef$mean
  MM<-Mdef$mean
     
  for(nn in 1:nsims){
     if(fmsymdef$dist!="none") fmsymX<-getrandom(1,fmsymdef$dist,a=fmsymdef$low,b=fmsymdef$up,mean=fmsymdef$mean,sd=fmsymdef$sd)
     if(btkdef$dist!="none") btkX<-getrandom(1, btkdef$dist,a= btkdef$low,b= btkdef$up,mean= btkdef$mean,sd= btkdef$sd)
     if(bmsykdef$dist!="none") bmsykX<-getrandom(1,bmsykdef$dist,a=bmsykdef$low,b=bmsykdef$up,mean=bmsykdef$mean,sd=bmsykdef$sd)
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
     f<-function(d){(bmsykX-(d^(1/(1-d))))^2}
      n<-optimize(f,c(0,1000),tol=0.000001)[[1]]
      g<-(n^(n/(n-1)))/(n-1)
     # Find K 
      findk<-function(K){ 
        B<-NULL;P<-NULL 
        Fmsy<-fmsymX*MM
        Umsy<-(Fmsy/(Fmsy+MM))*(1-exp(-Fmsy-MM))
        MSY<-K*bmsykX*Umsy
        B[1]<-K
       for(t in 1:timelen){  
         if(t<=agemat)  P[t]<-0
         if(t>agemat){
           if(bmsykX>=0.5) P[t]<-g*MSY*(B[t-agemat]/K)-g*MSY*(B[t-agemat]/K)^n 
           if(bmsykX>0.3 & bmsykX<0.5){
              bjoin<-(0.75*bmsykX-0.075)*K
              if(B[t-agemat]<bjoin){
               PJ<-g*MSY*(bjoin/K)-g*MSY*(bjoin/K)^n
               cc<-(1-n)*MSY*g*(bjoin^(n-2))*K^-n
               P[t]<-B[t-agemat]*(PJ/bjoin+cc*(B[t-agemat]-bjoin)) 
              }
              if(B[t-agemat]>=bjoin) P[t]<-g*MSY*(B[t-agemat]/K)-g*MSY*(B[t-agemat]/K)^n 
           }
           if(bmsykX<=0.3){
             bjoin<-(0.5*bmsykX)*K
              if(B[t-agemat]<bjoin){
               PJ<-g*MSY*(bjoin/K)-g*MSY*(bjoin/K)^n
               cc<-(1-n)*MSY*g*(bjoin^(n-2))*K^-n
               P[t]<-B[t-agemat]*(PJ/bjoin+cc*(B[t-agemat]-bjoin)) 
              }
              if(B[t-agemat]>=bjoin) P[t]<-g*MSY*(B[t-agemat]/K)-g*MSY*(B[t-agemat]/K)^n 
           }
         } 
        if(P[t]<0) P[t]<-0
         B[t+1]<-max(0,B[t]+P[t]-dcatch[t])
       }
       (btkX-(B[refyr]/K))^2 
  }
     
      out<-optimize(findk,c(kdef$low,kdef$up),tol=kdef$tol)
      bigK<-out$minimum
    
        B<-NULL;P<-NULL 
        Fmsy<-fmsymX*MM
        Umsy<-(Fmsy/(Fmsy+MM))*(1-exp(-Fmsy-MM))
        MSY<-bigK*bmsykX*Umsy
        B[1]<-bigK
       for(t in 1:timelen){       
         if(t<=agemat)  P[t]<-0
         if(t>agemat){
           if(bmsykX>=0.5) P[t]<-g*MSY*(B[t-agemat]/bigK)-g*MSY*(B[t-agemat]/bigK)^n 
           if(bmsykX>0.3 & bmsykX<0.5){
              bjoin<-(0.75*bmsykX-0.075)*bigK
              if(B[t-agemat]<bjoin){
               PJ<-g*MSY*(bjoin/bigK)-g*MSY*(bjoin/bigK)^n
               cc<-(1-n)*MSY*g*(bjoin^(n-2))*bigK^-n
               P[t]<-B[t-agemat]*(PJ/bjoin+cc*(B[t-agemat]-bjoin)) 
              }
              if(B[t-agemat]>=bjoin) P[t]<-g*MSY*(B[t-agemat]/bigK)-g*MSY*(B[t-agemat]/bigK)^n 
           }
           if(bmsykX<=0.3){
              bjoin<-(0.5*bmsykX)*bigK
              if(B[t-agemat]<bjoin){
                PJ<-g*MSY*(bjoin/bigK)-g*MSY*(bjoin/bigK)^n
                cc<-(1-n)*MSY*g*(bjoin^(n-2))*bigK^-n
                P[t]<-B[t-agemat]*(PJ/bjoin+cc*(B[t-agemat]-bjoin)) 
              }
                if(B[t-agemat]>=bjoin) P[t]<-g*MSY*(B[t-agemat]/bigK)-g*MSY*(B[t-agemat]/bigK)^n 
           }
         } 
         if(P[t]<0) P[t]<-0
         B[t+1]<-max(0,B[t]+P[t]-dcatch[t])
       }
       
      bll<-0
      if(min(B)>0 && max(B)<=bigK && (out$objective<=kdef$tol^2) && (abs((max(B)-bigK)/bigK)*100)<=kdef$permax) bll<-1 
      if(nn==1){
            write.table(t(c(bll,B)),file="Biotraj-dbsra.csv",sep=",",row.names=FALSE,col.names=FALSE,append=FALSE)
           if(catchout==1) write.table(t(c(bll,dcatch)),file="Catchtraj-dbsra.csv",sep=",",row.names=FALSE,col.names=FALSE,append=FALSE)
        }  
      if(nn>1){
        write.table(t(c(bll,B)),file="Biotraj-dbsra.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
        if(catchout==1) write.table(t(c(bll,dcatch)),file="Catchtraj-dbsra.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
      }
        storep[nn,1]<-bll
        storep[nn,2]<-fmsymX
        storep[nn,3]<-btkX
        storep[nn,4]<-bmsykX
        storep[nn,5]<-MM
        storep[nn,6]<-bigK
        storep[nn,7]<-Fmsy
        storep[nn,8]<-Umsy
        storep[nn,9]<-MSY
        storep[nn,10]<-bigK*bmsykX #Bmsy
        storep[nn,11]<-Umsy*B[timelen+1]  #OFL
        storep[nn,12]<-B[refyr]  
        storep[nn,13]<-B[timelen+1]
    }#n loop  
   
# Get results for bll==1 for graphing
  colnames(storep)<-c("ll","FmsyM","BtK","BmsyK","M","K","Fmsy","Umsy","MSY","Bmsy","OFLT1","Brefyr","BT1")
  datar<-storep[storep[,1]==1,]
if(length(datar[,1])>0){ 
    mMSY<-round(mean(datar$MSY),3)
    MSY95<-round(quantile(datar$MSY,probs=c(0.025,0.50,0.975)),4)
    mk<-round(mean(datar$K),3)
    k95<-round(quantile(datar$K,probs=c(0.025,0.50,0.975)),4)
    mBMSY<-round(mean(datar$Bmsy),3)
    BMSY95<-round(quantile(datar$Bmsy,probs=c(0.025,0.50,0.975)),4)
    mFMSY<-round(mean(datar$Fmsy),3)
    FMSY95<-round(quantile(datar$Fmsy,probs=c(0.025,0.50,0.975)),4)
    mUMSY<-round(mean(datar$Umsy),3)
    UMSY95<-round(quantile(datar$Umsy,probs=c(0.025,0.50,0.975)),4)
    mOFL<-round(mean(datar$OFLT1),3)
    OFL95<-round(quantile(datar$OFLT1,probs=c(0.025,0.50,0.975)),4)
    mM<-round(mean(datar$M),3)
    M95<-round(quantile(datar$M,probs=c(0.025,0.50,0.975)),4)
    mbtk<-round(mean(datar$BtK),3)
    btk95<-round(quantile(datar$BtK,probs=c(0.025,0.50,0.975)),4)
    mFmsyM<-round(mean(datar$FmsyM),3)
    FmsyM95<-round(quantile(datar$FmsyM,probs=c(0.025,0.50,0.975)),4)
    mBmsyK<-round(mean(datar$BmsyK),3)
    BmsyK95<-round(quantile(datar$BmsyK,probs=c(0.025,0.50,0.975)),4)
    mBrefyr<-round(mean(datar$Brefyr),3)
    Brefyr95<-round(quantile(datar$Brefyr,probs=c(0.025,0.50,0.975)),4)
   
    ### Graphs######################
   
    # Plot graphs
  if(grout>0){
     grunits<-data.frame(gr=c(1:14),cexa=0,cexl=0,cexx=0,nclass=0,mains=" ",cexmain=0,lwd=0,stringsAsFactors=FALSE)
     grunits$lwd<-ifelse(grunits$gr %in% c(1,13),grdef$lwd,0)   
     grunits[,2]<-grdef$cex.axis
     grunits[,3]<-grdef$cex.lab
     grunits[,4]<-grdef$cex
     grunits$nclass<-ifelse(grunits$gr %in% c(2,3,4,5,6,7,8,9,10,11,12,14),grdef$nclasses,grunits$nclass)
     grunits$cexmain<-grdef$cex.main
     grunits[graphs,6]<-grdef$mains
    
    if(any(graphs==1)){
          plot(catch~year,type="l",xlab="Year",ylab=paste("Catch (",catdef$unit,")",sep=""),
          ylim=c(0,round(max(catch,mMSY,MSY95[3]),0)),cex=grunits[1,4],cex.lab=grunits[1,3],
            cex.axis=grunits[1,2],main=grunits[1,6],cex.main=grunits[1,7],lwd=grunits[1,8])
           if(pstdef$ol==1){
             abline(h=MSY95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(h=MSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
             abline(h=MSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
           }
      
          if(grout==2){
           word.tif("catch")
            plot(catch~year,type="l",xlab="Year",ylab=paste("Catch (",catdef$unit,")",sep=""),
            ylim=c(0,round(max(catch,mMSY,MSY95[3]),0)),cex=grunits[1,4],cex.lab=grunits[1,3],
              cex.axis=grunits[1,2],main=grunits[1,6],cex.main=grunits[1,7],lwd=grunits[1,8])
             if(pstdef$ol==1){
             abline(h=MSY95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(h=MSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
             abline(h=MSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
             }
           dev.off()
         }
      }
      if(any(graphs==2)){
        hist(datar$K,freq=FALSE,xlim=c(0,max(datar$K)*1.40),xlab=paste("K (",catdef$unit,")",sep=""),nclass=grunits[3,5],cex.lab=grunits[3,3],cex.axis=grunits[3,2],main=grunits[3,6],cex.main=grunits[3,7])
           if(pstdef$ol==1){ 
             abline(v=k95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(v=k95[1],lwd=pstdef$lwd,lty=pstdef$llty)
             abline(v=k95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
            }
        if(grout==2){
            word.tif("Kden") 
             hist(datar$K,freq=FALSE,xlim=c(0,max(datar$K)*1.40),xlab=paste("K (",catdef$unit,")",sep=""),nclass=grunits[3,5],cex.lab=grunits[3,3],cex.axis=grunits[3,2],main=grunits[3,6],cex.main=grunits[3,7])
           if(pstdef$ol==1){ 
             abline(v=k95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(v=k95[1],lwd=pstdef$lwd,lty=pstdef$llty)
             abline(v=k95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
            }
          dev.off()
         }
      }
      if(any(graphs==3)){
           hist(datar$Bmsy,freq=FALSE,xlim=c(0,max(datar$Bmsy)*1.40),xlab=paste("Bmsy (",catdef$unit,")",sep=""),nclass=grunits[4,5],cex.lab=grunits[4,3],cex.axis=grunits[4,2],main=grunits[4,6],cex.main=grunits[4,7])
            if(pstdef$ol==1){ 
             abline(v=BMSY95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(v=BMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
             abline(v=BMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
            }
        if(grout==2){
           word.tif("Bmsyden") 
           hist(datar$Bmsy,freq=FALSE,xlim=c(0,max(datar$Bmsy)*1.40),xlab=paste("Bmsy (",catdef$unit,")",sep=""),nclass=grunits[4,5],cex.lab=grunits[4,3],cex.axis=grunits[4,2],main=grunits[4,6],cex.main=grunits[4,7])
            if(pstdef$ol==1){ 
             abline(v=BMSY95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
             abline(v=BMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
             abline(v=BMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
            }
          dev.off()
         }
      }    
      if(any(graphs==4)){      
          hist(datar$MSY,freq=FALSE,xlim=c(min(datar$MSY)*0.80,max(datar$MSY)*1.20),xlab=paste("MSY (",catdef$unit,")",sep=""),nclass=grunits[5,5],cex.lab=grunits[5,3],cex.axis=grunits[5,2],main=grunits[5,6],cex.main=grunits[5,7])
          if(pstdef$ol==1){ 
           abline(v=MSY95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
           abline(v=MSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
           abline(v=MSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
          }
         if(grout==2){
           word.tif("MSYden") 
            hist(datar$MSY,freq=FALSE,xlim=c(min(datar$MSY)*0.80,max(datar$MSY)*1.20),xlab=paste("MSY (",catdef$unit,")",sep=""),nclass=grunits[5,5],cex.lab=grunits[5,3],cex.axis=grunits[5,2],main=grunits[5,6],cex.main=grunits[5,7])
          if(pstdef$ol==1){ 
           abline(v=MSY95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
           abline(v=MSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
           abline(v=MSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
          }
           dev.off()
          }
      }  
    if(any(graphs==5)){      
          hist(datar$Fmsy,freq=FALSE,xlim=c(0,max(datar$Fmsy)*1.30),xlab="Fmsy",nclass=grunits[6,5],cex.lab=grunits[6,3],cex.axis=grunits[6,2],main=grunits[6,6],cex.main=grunits[6,7])
          if(pstdef$ol==1){ 
           abline(v=FMSY95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
           abline(v=FMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
           abline(v=FMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
          }
         if(grout==2){
           word.tif("Fmsyden") 
            hist(datar$Fmsy,freq=FALSE,xlim=c(0,max(datar$Fmsy)*1.20),xlab="Fmsy",nclass=grunits[6,5],cex.lab=grunits[6,3],cex.axis=grunits[6,2],main=grunits[6,6],cex.main=grunits[6,7])
          if(pstdef$ol==1){ 
           abline(v=FMSY95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
           abline(v=FMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
           abline(v=FMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
          }
           dev.off()
          }
      }  
     if(any(graphs==6)){      
         hist(datar$Umsy,freq=FALSE,xlim=c(0,max(datar$Umsy)*1.20),xlab="Umsy",nclass=grunits[7,5],cex.lab=grunits[7,3],cex.axis=grunits[7,2],main=grunits[7,6],cex.main=grunits[7,7])
         if(pstdef$ol==1){ 
           abline(v=UMSY95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
           abline(v=UMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
           abline(v=UMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
        if(grout==2){
           word.tif("Umsyden") 
              hist(datar$Umsy,freq=FALSE,xlim=c(0,max(datar$Umsy)*1.20),xlab="Umsy",nclass=grunits[7,5],cex.lab=grunits[7,3],cex.axis=grunits[7,2],main=grunits[7,6],cex.main=grunits[7,7])
            if(pstdef$ol==1){ 
            abline(v=UMSY95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
            abline(v=UMSY95[1],lwd=pstdef$llwd,lty=pstdef$llty)
            abline(v=UMSY95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
            }
           dev.off()
          }
     }
    if(any(graphs==7)){   
       hist(datar$OFL,freq=FALSE,xlim=c(min(datar$OFL)*0.8,max(datar$OFL)*1.30),xlab=paste("OFL (",catdef$unit,")",sep=""),nclass=grunits[8,5],cex.lab=grunits[8,3],cex.axis=grunits[8,2],main=grunits[8,6],cex.main=grunits[8,7])
        if(pstdef$ol==1){
          abline(v=OFL95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=OFL95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=OFL95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
        }
      if(grout==2){
           word.tif("OFLden") 
           hist(datar$OFL,freq=FALSE,xlim=c(min(datar$OFL)*0.8,max(datar$OFL)*1.30),xlab=paste("OFL (",catdef$unit,")",sep=""),nclass=grunits[8,5],cex.lab=grunits[8,3],cex.axis=grunits[8,2],main=grunits[8,6],cex.main=grunits[8,7])
           if(pstdef$ol==1){
            abline(v=OFL95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
            abline(v=OFL95[1],lwd=pstdef$llwd,lty=pstdef$llty)
            abline(v=OFL95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
           }
           dev.off()
       }
    }  
    if(any(graphs==8)){   
         hist(datar$M,freq=FALSE,xlim=c(0,max(datar$M)*1.20),xlab="M",nclass=grunits[9,5],cex.lab=grunits[9,3],cex.axis=grunits[9,2],main=grunits[9,6],cex.main=grunits[9,7])
         if(pstdef$ol==1){
          abline(v=M95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=M95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=M95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
        if(grout==2){
           word.tif("Mden") 
           hist(datar$M,freq=FALSE,xlim=c(0,max(datar$M)*1.20),xlab="M",nclass=grunits[9,5],cex.lab=grunits[9,3],cex.axis=grunits[9,2],main=grunits[9,6],cex.main=grunits[9,7])
         if(pstdef$ol==1){
          abline(v=M95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=M95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=M95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
          dev.off()
        }
    }
    if(any(graphs==9)){   
         hist(datar$BtK,freq=FALSE,xlim=c(0,max(datar$BtK)*1.30),xlab="Bt/K",nclass=grunits[10,5],cex.lab=grunits[10,3],cex.axis=grunits[10,2],main=grunits[10,6],cex.main=grunits[10,7])
         if(pstdef$ol==1){
          abline(v=btk95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=btk95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=btk95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
        if(grout==2){
           word.tif("BtKden") 
           hist(datar$BtK,freq=FALSE,xlim=c(0,max(datar$BtK)*1.30),xlab="Bt/K",nclass=grunits[10,5],cex.lab=grunits[10,3],cex.axis=grunits[10,2],main=grunits[10,6],cex.main=grunits[10,7])
         if(pstdef$ol==1){
          abline(v=btk95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=btk95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=btk95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
          dev.off()
        }
    }
       if(any(graphs==10)){   
         hist(datar$FmsyM,freq=FALSE,xlim=c(0,max(datar$FmsyM)*1.30),xlab="Fmsy/M",nclass=grunits[11,5],cex.lab=grunits[11,3],cex.axis=grunits[11,2],main=grunits[11,6],cex.main=grunits[11,7])
         if(pstdef$ol==1){
          abline(v=FmsyM95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=FmsyM95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=FmsyM95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
        if(grout==2){
           word.tif("FmsyMden") 
           hist(datar$FmsyM,freq=FALSE,xlim=c(0,max(datar$FmsyM)*1.30),xlab="Fmsy/M",nclass=grunits[11,5],cex.lab=grunits[11,3],cex.axis=grunits[11,2],main=grunits[11,6],cex.main=grunits[11,7])
         if(pstdef$ol==1){
          abline(v=FmsyM95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=FmsyM95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=FmsyM95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
          dev.off()
        }
    }
     if(any(graphs==11)){   
         hist(datar$BmsyK,freq=FALSE,xlim=c(0,max(datar$BmsyK)*1.30),xlab="Bmsy/K",nclass=grunits[11,5],cex.lab=grunits[11,3],cex.axis=grunits[11,2],main=grunits[11,6],cex.main=grunits[11,7])
         if(pstdef$ol==1){
          abline(v=BmsyK95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=BmsyK95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=BmsyK95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
        if(grout==2){
           word.tif("BmsyKden") 
          hist(datar$BmsyK,freq=FALSE,xlim=c(0,max(datar$BmsyK)*1.30),xlab="Bmsy/K",nclass=grunits[11,5],cex.lab=grunits[11,3],cex.axis=grunits[11,2],main=grunits[11,6],cex.main=grunits[11,7])
         if(pstdef$ol==1){
          abline(v=BmsyK95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=BmsyK95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=BmsyK95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
          dev.off()
        }
    }
     if(any(graphs==12)){   
         hist(datar$Brefyr,freq=FALSE,xlim=c(min(datar$Brefyr)*0.7,max(datar$Brefyr)*1.30),xlab="Brefyr",nclass=grunits[11,5],cex.lab=grunits[11,3],cex.axis=grunits[11,2],main=grunits[11,6],cex.main=grunits[11,7])
         if(pstdef$ol==1){
          abline(v=Brefyr95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=Brefyr95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=Brefyr95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
        if(grout==2){
           word.tif("Brefyrden") 
          hist(datar$Brefyr,freq=FALSE,xlim=c(min(datar$Brefyr)*0.7,max(datar$Brefyr)*1.30),xlab="Brefyr",nclass=grunits[11,5],cex.lab=grunits[11,3],cex.axis=grunits[11,2],main=grunits[11,6],cex.main=grunits[11,7])
         if(pstdef$ol==1){
          abline(v=Brefyr95[2],lwd=pstdef$mlwd,lty=pstdef$mlty)
          abline(v=Brefyr95[1],lwd=pstdef$llwd,lty=pstdef$llty)
          abline(v=Brefyr95[3],lwd=pstdef$ulwd,lty=pstdef$ulty)
         }
          dev.off()
        }
    }
     if(any(graphs==13)){
         bigB<-read.csv("Biotraj-dbsra.csv",header=FALSE)
         ar<-bigB[,1]
         bigB<-t(bigB[,-1])
         cols<-ifelse(ar==0,"gray85","black")
         types<-ifelse(ar==0,3,1)
                 
          par(mfrow=c(1,2))
         matplot(y=bigB[,c(which(ar==1))],x=c(year,year[length(year)]+1),type="l",lty=types[c(which(ar==1))],xlab="Year",ylab=paste("Biomass (",catdef$unit,")",sep=""),
            ylim=c(0,max(bigB)),cex=grunits[13,4],cex.lab=grunits[13,3],col=cols[c(which(ar==1))],main="Accepted",
              cex.axis=grunits[13,2],cex.main=grunits[13,7],lwd=grunits[13,8])    
           lines(y=apply(bigB[,c(which(ar==1))],1,median),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$mlty,col="red")
           lines(y=apply(bigB[,c(which(ar==1))],1,function(x){quantile(x,probs=0.025)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$llty,col="red")
           lines(y=apply(bigB[,c(which(ar==1))],1,function(x){quantile(x,probs=0.975)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$ulty,col="red")
 
         matplot(y=bigB[,c(which(ar==0))],x=c(year,year[length(year)]+1),type="l",lty=types[c(which(ar==0))],xlab="Year",ylab=paste("Biomass (",catdef$unit,")",sep=""),
            ylim=c(0,max(bigB)),cex=grunits[13,4],cex.lab=grunits[13,3],col=cols[c(which(ar==0))],main="Rejected",
              cex.axis=grunits[13,2],cex.main=grunits[13,7],lwd=grunits[13,8])    
           lines(y=apply(bigB[,c(which(ar==0))],1,median),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$mlty,col="red")
           lines(y=apply(bigB[,c(which(ar==0))],1,function(x){quantile(x,probs=0.025)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$llty,col="red")
           lines(y=apply(bigB[,c(which(ar==0))],1,function(x){quantile(x,probs=0.975)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$ulty,col="red")
        if(grout==2){
            word.tif("Biomasstraj") 
            par(mfrow=c(1,2))
            matplot(y=bigB[,c(which(ar==1))],x=c(year,year[length(year)]+1),type="l",lty=types[c(which(ar==1))],xlab="Year",ylab=paste("Biomass (",catdef$unit,")",sep=""),
            ylim=c(0,max(bigB)),cex=grunits[13,4],cex.lab=grunits[13,3],col=cols[c(which(ar==1))],main="Accepted",
              cex.axis=grunits[13,2],cex.main=grunits[13,7],lwd=grunits[13,8])    
           lines(y=apply(bigB[,c(which(ar==1))],1,median),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$mlty,col="red")
           lines(y=apply(bigB[,c(which(ar==1))],1,function(x){quantile(x,probs=0.025)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$llty,col="red")
           lines(y=apply(bigB[,c(which(ar==1))],1,function(x){quantile(x,probs=0.975)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$ulty,col="red")
 
         matplot(y=bigB[,c(which(ar==0))],x=c(year,year[length(year)]+1),type="l",lty=types[c(which(ar==0))],xlab="Year",ylab=paste("Biomass (",catdef$unit,")",sep=""),
            ylim=c(0,max(bigB)),cex=grunits[13,4],cex.lab=grunits[13,3],col=cols[c(which(ar==0))],main="Rejected",
              cex.axis=grunits[13,2],cex.main=grunits[13,7],lwd=grunits[13,8])    
           lines(y=apply(bigB[,c(which(ar==0))],1,median),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$mlty,col="red")
           lines(y=apply(bigB[,c(which(ar==0))],1,function(x){quantile(x,probs=0.025)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$llty,col="red")
           lines(y=apply(bigB[,c(which(ar==0))],1,function(x){quantile(x,probs=0.975)}),x=c(year,year[length(year)]+1),lwd=2,lty=pstdef$ulty,col="red")
          dev.off()
        }
       rm(bigB)
       par(mfrow=c(1,1))
    }
         if(any(graphs==14)){ 
          #BmsyK
           breaks<-hist(storep[,"BmsyK"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"BmsyK"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"BmsyK"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="Bmsy/K",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
         #M
           breaks<-hist(storep[,"M"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"M"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"M"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="M",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
         #BtK
           breaks<-hist(storep[,"BtK"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"BtK"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"BtK"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="Bt/K",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
            #FmsyM
           breaks<-hist(storep[,"FmsyM"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"FmsyM"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"FmsyM"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="Fmsy/M",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE)  
           
           #K
           breaks<-hist(storep[,"K"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"K"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"K"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="K",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
          
            #MSY
           breaks<-hist(storep[,"MSY"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"MSY"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"MSY"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="MSY",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
           
          #Bmsy
           breaks<-hist(storep[,"Bmsy"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"Bmsy"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"Bmsy"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="Bmsy",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
           
           #Umsy
           breaks<-hist(storep[,"Umsy"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"Umsy"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"Umsy"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="Umsy",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
           if(grout==2){
           word.tif("BmsykAR") 
           breaks<-hist(storep[,"BmsyK"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"BmsyK"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"BmsyK"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="Bmsy/K",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
          dev.off()
           word.tif("MAR") 
             breaks<-hist(storep[,"M"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"M"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"M"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="M",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
            dev.off()
           word.tif("BtKAR")
              #BtK
           breaks<-hist(storep[,"BtK"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"BtK"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"BtK"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="Bt/K",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
      dev.off()
          word.tif("FmsyMAR")
           #FmsyM
           breaks<-hist(storep[,"FmsyM"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"FmsyM"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"FmsyM"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="Fmsy/M",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
             dev.off()
            #K
            word.tif("KAR")
              breaks<-hist(storep[,"K"],plot=FALSE,nclass=grunits[14,5])$breaks
             good<-table((cut(storep[storep[,1]==1,"K"], breaks)))
             dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
             bad<-table((cut(storep[storep[,1]==0,"K"], breaks)))
             dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
             newd<-cbind(good,bad)       
             xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
             barplot(t(newd),xlab="K",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
             axis(1,at=xpos,labels=FALSE)
            dev.off()
             #MSY
             
            #MSY
           word.tif("MSYAR")
           breaks<-hist(storep[,"MSY"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"MSY"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"MSY"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="MSY",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
           dev.off()
             #Bmsy
          word.tif("BmsyAR")
           breaks<-hist(storep[,"Bmsy"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"Bmsy"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"Bmsy"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="Bmsy",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE)
          dev.off()
        
          word.tif("UmsyAR")
           breaks<-hist(storep[,"Umsy"],plot=FALSE,nclass=grunits[14,5])$breaks
           good<-table((cut(storep[storep[,1]==1,"Umsy"], breaks)))
           dimnames(good)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           bad<-table((cut(storep[storep[,1]==0,"Umsy"], breaks)))
           dimnames(bad)[[1]]<-(breaks[1:c(length(breaks)-1)]+breaks[2:c(length(breaks))])/2    
           newd<-cbind(good,bad)       
           xpos<-barplot(t(newd),ylim=c(0,max(newd)),plot=FALSE)  
           barplot(t(newd),xlab="Umsy",ylab="Frequency",col=c("black","NA"),cex.main=0.7,main=paste("Accepted=Black","  ","Rejected=White"))
           axis(1,at=xpos,labels=FALSE) 
          dev.off()
        }
    }
  }#grout>0

  ### Output
    outs<-data.frame(Distr=NA,Min=NA,Max=NA,Mean=NA,sd=NA)
    outs[1,]<-cbind(fmsymdef$dist,fmsymdef$low,fmsymdef$up,fmsymdef$mean,fmsymdef$sd)
    outs[2,]<-cbind(btkdef$dist,btkdef$low,btkdef$up,btkdef$mean,btkdef$sd)
    outs[3,]<-cbind(bmsykdef$dist,bmsykdef$low,bmsykdef$up,bmsykdef$mean,bmsykdef$sd)
    outs[4,]<-cbind(Mdef$dist,Mdef$low,Mdef$up,Mdef$mean,Mdef$sd)
    outs[5,]<-cbind(NA,btkdef$refyr,NA,NA,NA)
    colnames(outs)<-c("Distr","Lower","Upper","Mean ","SD")
    rownames(outs)<-c("Fmsy/M","Br/K","Bmsy/K","M","refyr")
  
    outs1<-data.frame(Mean=NA,Median=NA,per2_5=NA,per97_5=NA,min=NA,max=NA)
    outs1[1,]<-cbind(mFmsyM,FmsyM95[[2]],FmsyM95[[1]],FmsyM95[[3]],min(datar$FmsyM),max(datar$FmsyM))
    outs1[2,]<-cbind(mbtk,btk95[[2]],btk95[[1]],btk95[[3]],min(datar$BtK),max(datar$BtK))
    outs1[3,]<-cbind(mBmsyK,BmsyK95[[2]],BmsyK95[[1]],BmsyK95[[3]],min(datar$BmsyK),max(datar$BmsyK))
    outs1[4,]<-cbind(mM,M95[[2]],M95[[1]],M95[[3]],min(datar$M),max(datar$M))
    colnames(outs1)<-c("Mean (ll=1)","Median (ll=1)","2.5% (ll=1)","97.5% (ll=1)","min (ll=1)","max (ll=1)")
    rownames(outs1)<-c("Fmsy/M","Bt/K","Bmsy/K","M")
    
    outs2<-data.frame(Mean=NA,Median=NA, per2_5=NA,per97_5=NA,min=NA,max=NA)
    outs2[1,]<-cbind(mMSY,MSY95[[2]],MSY95[[1]],MSY95[[3]],min(datar$MSY),max(datar$MSY))
    outs2[2,]<-cbind(mBMSY,BMSY95[[2]],BMSY95[[1]],BMSY95[[3]],min(datar$Bmsy),max(datar$Bmsy))
    outs2[3,]<-cbind(mFMSY,FMSY95[[2]],FMSY95[[1]],FMSY95[[3]],min(datar$Fmsy),max(datar$Fmsy))
    outs2[4,]<-cbind(mUMSY,UMSY95[[2]],UMSY95[[1]],UMSY95[[3]],min(datar$Umsy),max(datar$Umsy))
    outs2[5,]<-cbind(mOFL,OFL95[[2]],OFL95[[1]],OFL95[[3]],min(datar$OFL),max(datar$OFL))
    outs2[6,]<-cbind(mBrefyr,Brefyr95[[2]],Brefyr95[[1]],Brefyr95[[3]],min(datar$Brefyr),max(datar$Brefyr))
    outs2[7,]<-cbind(mk,k95[[2]],k95[[1]],k95[[3]],min(datar$K),max(datar$K))
    colnames(outs2)<-c("Mean (ll=1)","Median (ll=1)","2.5% (ll=1)","97.5% (ll=1)","min (ll=1)","max (ll=1)")
    rownames(outs2)<-c("MSY","Bmsy","Fmsy","Umsy","OFL","Brefyr","K")
    
      ans<-list(outs,outs1,outs2,storep,agemat,c(max(year)+1),"dbsra")
      names(ans)<-c("Initial","Parameters","Estimates","Values","agemat","end1yr","type")
    
  } #length(datar[,1])>1
  
if(length(datar[,1])==0){
       outs<-data.frame(Distr=NA,Min=NA,Max=NA,Mean=NA,sd=NA)
    outs[1,]<-cbind(fmsymdef$dist,fmsymdef$low,fmsymdef$up,fmsymdef$mean,fmsymdef$sd)
    outs[2,]<-cbind(btkdef$dist,btkdef$low,btkdef$up,btkdef$mean,btkdef$sd)
    outs[3,]<-cbind(bmsykdef$dist,bmsykdef$low,bmsykdef$up,bmsykdef$mean,bmsykdef$sd)
    outs[4,]<-cbind(Mdef$dist,Mdef$low,Mdef$up,Mdef$mean,Mdef$sd)
    outs[5,]<-cbind(NA,btkdef$refyr,NA,NA,NA)
    colnames(outs)<-c("Distr","Lower","Upper","Mean ","SD")
    rownames(outs)<-c("Fmsy/M","Bt/K","Bmsy/K","M","refyr")
    outs1<-data.frame(Mean=NA,Median=NA,per2_5=NA,per97_5=NA,min=NA,max=NA)
    colnames(outs1)<-c("Mean (ll=1)","Median (ll=1)","2.5% (ll=1)","97.5% (ll=1)")
    outs2<-data.frame(Mean=NA,Median=NA, per2_5=NA,per97_5=NA,min=NA,max=NA)
    colnames(outs2)<-c("Mean (ll=1)","Median (ll=1)","2.5% (ll=1)","97.5% (ll=1)","min (ll=1)","max (ll=1)")
      ans<-list(outs,outs1,outs2,storep,agemat,c(max(year)+1),"dbsra")
      names(ans)<-c("Initial","Parameters","Estimates","Values","agemat","end1yr","type")
   warning("None of the runs had a likelihood equal to 1")
  }
    return(ans) 
}#function

