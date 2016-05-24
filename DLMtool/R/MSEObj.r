
# A collection of functions which operate on MSE object (class MSE) which is 
# produced by the DLMtool MSE.
# These functions are used to summarize, plot, and manipulate the results of the
# MSE. 
# Most of these functions have accompanying help files (?function.name), 
# although, for purpose of clarity, some supporting internal functions have also
# been included here.

# January 2016
# Tom Carruthers UBC (t.carruthers@fisheries.ubc.ca)
# Adrian Hordyk (a.hordyk@murdoch.edu.au)

# Have I undertaken enough simulations (nsim)? 
# Has my MSE converged on stable (reliable) peformance metrics?
CheckConverg <- function(MSEobj, thresh=2, Plot=TRUE) {
  nm<-MSEobj@nMPs
  nsim<-MSEobj@nsim
  proyears<-MSEobj@proyears
  
  Yd <- CumlYd <- array(NA,c(nm,nsim))
  P10 <- CumlP10 <- array(NA,c(nm,nsim))
  P50 <- CumlP50 <- array(NA,c(nm,nsim))
  P100 <- CumlP100 <- array(NA,c(nm,nsim))
  POF <- CumlPOF <- array(NA,c(nm,nsim))
  yind <-max(MSEobj@proyears-4,1):MSEobj@proyears
  RefYd<-MSEobj@OM$RefY
  
  for(m in 1:nm){
    Yd[m,] <-round(apply(MSEobj@C[,m,yind],1,mean,na.rm=T)/RefYd*100,1)
    POF[m,] <-round(apply(MSEobj@F_FMSY[,m,]>1,1,sum,na.rm=T)/proyears*100,1)
    P10[m,] <-round(apply(MSEobj@B_BMSY[,m,]<0.1,1,sum,na.rm=T)/proyears*100,1)
    P50[m,] <-round(apply(MSEobj@B_BMSY[,m,]<0.5,1,sum,na.rm=T)/proyears*100,1)
    P100[m,] <-round(apply(MSEobj@B_BMSY[,m,]<1,1,sum,na.rm=T)/proyears*100,1)
    CumlYd[m,] <- cumsum(Yd[m,]) / seq_along(Yd[m,])#/ mean(Yd[m,], na.rm=TRUE) 
    CumlPOF[m,] <- cumsum(POF[m,]) / seq_along(POF[m,])# / mean(POF[m,], na.rm=TRUE)
    CumlP10[m,] <- cumsum(P10[m,]) / seq_along(P10[m,])# / mean(P10[m,], na.rm=TRUE)
    CumlP50[m,] <- cumsum(P50[m,]) / seq_along(P50[m,])# / mean(P50[m,], na.rm=TRUE)
    CumlP100[m,] <- cumsum(P100[m,]) / seq_along(P100[m,])# / mean(P100[m,], na.rm=TRUE)
  }
  
  # CumlYd[is.nan(CumlYd)] <- 1
  # CumlPOF[is.nan(CumlPOF)] <- 1
  # CumlP10[is.nan(CumlP10)] <- 1
  # CumlP50[is.nan(CumlP50)] <- 1
  # CumlP100[is.nan(CumlP100)] <- 1
  if (Plot) {
    par(mfrow=c(2,3), cex.axis=1.5, cex.lab=1.7, oma=c(1,1,0,0), mar=c(5,5,1,1), bty="l")
    matplot(t(CumlYd), type="l", xlab="Iteration", ylab="Rel. Yield")
    matplot(t(CumlPOF), type="l", xlab="Iteration", ylab="Prob. F/FMSY > 1")
    matplot(t(CumlP10), type="l", xlab="Iteration", ylab="Prob. B/BMSY < 0.1")
    matplot(t(CumlP50), type="l", xlab="Iteration", ylab="Prob. B/BMSY < 0.5")
    matplot(t(CumlP100), type="l", xlab="Iteration", ylab="Prob. B/BMSY < 1")
  }
  
  Chk <- function(X) {
    # checks if difference in
    # last 10 iterations is greater than thresh
    L <- length(X)
    Y <- 1: min(nsim, 10)
    
    # return(mean(abs((X[L-1:10] - X[L]))/X[L], na.rm=TRUE) > thresh)
    return(mean(abs((X[L-Y] - X[L])), na.rm=TRUE) > thresh)
  }
  
  NonCon <- sort(unique(c(which(apply(CumlYd, 1, Chk)), which(apply(CumlPOF, 1, Chk)),
                          which(apply(CumlP10, 1, Chk)), which(apply(CumlP50, 1, Chk)), 	
                          which(apply(CumlP100, 1, Chk)))))
  
  if (length(NonCon) == 1) NonCon <- rep(NonCon, 2)
  if (length(NonCon) > 0) {
    if (Plot) {
      plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
      text(0.5, 0.5, "Some MPs have not converged", cex=1)
      #ev.new()
      par(mfrow=c(2,3), cex.axis=1.5, cex.lab=1.7, oma=c(1,1,0,0), mar=c(5,5,1,1), bty="l")
      matplot(t(CumlYd[NonCon,]), type="b", xlab="Iteration", ylab="Rel. Yield", lwd=2)
      matplot(t(CumlPOF[NonCon,]), type="b", xlab="Iteration", ylab="Prob. F/FMSY > 1", lwd=2)
      matplot(t(CumlP10[NonCon,]), type="b", xlab="Iteration", ylab="Prob. B/BMSY < 0.1", lwd=2)
      matplot(t(CumlP50[NonCon,]), type="b", xlab="Iteration", ylab="Prob. B/BMSY < 0.5", lwd=2)
      matplot(t(CumlP100[NonCon,]), type="b", xlab="Iteration", ylab="Prob. B/BMSY < 1", lwd=2)
      legend(nsim*1.25, 50, legend=MSEobj@MPs[NonCon], col=1:length(NonCon), bty="n", 
             xpd=NA, lty=1, lwd=2, cex=1.25)
    }  
    
    message("Some MPs may not have converged in ", nsim, " iterations (threshold = ", 
            thresh, "%)")
    message("MPs are: ", paste(MSEobj@MPs[NonCon], " "))
    message("MPs #: ", paste(NonCon, " "))
    return(data.frame(Num=NonCon, MP=MSEobj@MPs[NonCon]))
  }
  if (length(NonCon) == 0) {
    if (Plot) {
      plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
      text(0.5, 0.5, "All MPs converged", cex=1)
    }	
    message("All MPs appear to have converged in ", nsim, " iterations (threshold = ", 
            thresh, "%)")
  }

}

# Trade-Off Plots --------------------------------------------------------------
Tplot<-function(MSEobj,nam=NA){
  FMSYr<-quantile(MSEobj@F_FMSY,c(0.001,0.90),na.rm=T)
  BMSYr<-quantile(MSEobj@B_BMSY,c(0.001,0.975),na.rm=T)

  colsse<-rainbow(100,start=0,end=0.36)[1:100]
  colB<-rep(colsse[100],ceiling(BMSYr[2]*100))
  colB[1:100]<-colsse
  colB<-makeTransparent(colB,60)
  colsse<-rainbow(200,start=0,end=0.36)[200:1]
  colF<-rep(colsse[200],ceiling(FMSYr[2]*100))
  colF[1:200]<-colsse
  colF<-makeTransparent(colF,60)

  Yd<-rep(NA,MSEobj@nMPs)
  P10<-rep(NA,MSEobj@nMPs)
  P50<-rep(NA,MSEobj@nMPs)
  P100<-rep(NA,MSEobj@nMPs)
  POF<-rep(NA,MSEobj@nMPs)
  yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  RefYd<-MSEobj@OM$RefY

  for(mm in 1:MSEobj@nMPs){
    Yd[mm]<-round(mean(apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd,na.rm=T)*100,1)
  #cbind(MSEobj@C[,mm,yind],unlist(MSEobj@OM$MSY))
    POF[mm]<-round(sum(MSEobj@F_FMSY[,mm,]>1,na.rm=T)/prod(dim(MSEobj@F_FMSY[,mm,]),na.rm=T)*100,1)
    P10[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    P50[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    P100[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
  }
  
  #dev.new2(width=7,height=7)
  par(mfrow=c(2,2),mai=c(0.9,1,0.1,0.1),omi=c(0.1,0.1,0.4,0))

  tradeoffplot(POF,Yd,"Prob. of overfishing (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
  tradeoffplot(P100,Yd,"Prob. biomass < BMSY (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
  tradeoffplot(P50,Yd,"Prob. biomass < 0.5BMSY (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
  tradeoffplot(P10,Yd,"Prob. biomass < 0.1BMSY (%)", "Relative yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)

  if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  if(!is.na(nam) & !is.character(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
  if(!is.na(nam) & is.character(nam))mtext(nam,3,outer=T,line=0.3,font=2)
}

Tplot2<-function(MSEobj,nam=NA){
  LTY<-rep(NA,MSEobj@nMPs)
  STY<-rep(NA,MSEobj@nMPs)
  VY<-rep(NA,MSEobj@nMPs)
  B10<-rep(NA,MSEobj@nMPs)
  yend<-max(MSEobj@proyears-4,1):MSEobj@proyears
  ystart<-1:5
  RefYd<-MSEobj@OM$RefY
  y1<-1:(MSEobj@proyears-1)
  y2<-2:MSEobj@proyears
  for(mm in 1:MSEobj@nMPs){
    LTY[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)*100
    STY[mm]<-round(sum(MSEobj@C[,mm,ystart]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(ystart)),3)*100
    AAVY<-apply(((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T)/apply(MSEobj@C[,mm,y2],1,mean,na.rm=T)
    VY[mm]<-round(sum(AAVY<0.1,na.rm=T)/MSEobj@nsim,3)*100
    B10[mm]<-round(sum(MSEobj@B_BMSY[,mm,]>0.1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,])),3)*100
  }
  par(mfrow=c(1,2),mai=c(1.5,1.5,0.1,0.1),omi=c(0.1,0.1,0.4,0))
  tradeoffplot(STY,LTY,"P(Short term yield over half FMSY)","P(Long term yield over half FMSY)",MSEobj@MPs[1:MSEobj@nMPs],vl=1,hl=1)
  tradeoffplot(B10,VY,"P(Biomass > 0.1 BMSY)", "P(CV in yield less than 0.1)",MSEobj@MPs[1:MSEobj@nMPs],vl=1,hl=1)
  if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  if(!is.na(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
}

NOAA_plot<-function(MSEobj,nam=NA,type=NA,panel=T){
  
  Yd<-rep(NA,MSEobj@nMPs)
  B50<-rep(NA,MSEobj@nMPs)
  PNOF<-rep(NA,MSEobj@nMPs)
  LTY<-rep(NA,MSEobj@nMPs)
  STY<-rep(NA,MSEobj@nMPs)
  VY<-rep(NA,MSEobj@nMPs)
  
  y1<-1:(MSEobj@proyears-1)
  y2<-2:MSEobj@proyears
  
  yend<-max(MSEobj@proyears-4,1):MSEobj@proyears
  
  RefYd<-MSEobj@OM$RefY
  
  for(mm in 1:MSEobj@nMPs){
    
    PNOF[mm]<-round(sum(MSEobj@F_FMSY[,mm,]<1,na.rm=T)/prod(dim(MSEobj@F_FMSY[,mm,]),na.rm=T)*100,1)
    B50[mm]<-round(sum(MSEobj@B_BMSY[,mm,]>0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    LTY[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)*100
    AAVY<-apply((((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])/MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T) 
    VY[mm]<-round(sum(AAVY<0.15,na.rm=T)/MSEobj@nsim,3)*100

  }
  
  #dev.new2(width=7,height=7)
  if(panel)par(mfrow=c(1,2),mai=c(1.5,1.5,0.1,0.1),omi=c(0.1,0.1,0.4,0))
  
  if(is.na(type)){
    tradeoffplot(PNOF,LTY,"Prob. of not overfishing (%)", "Long-term yield ",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100)
    tradeoffplot(B50,VY,"Prob. biomass above half BMSY (%)", "Prob. AAVY less than 15%",MSEobj@MPs[1:MSEobj@nMPs],vl=80,hl=50)
  }else{
    tradeoffplot3(PNOF,LTY,"Prob. of not overfishing (%)", "Long-term yield",MSEobj@MPs[1:MSEobj@nMPs],vl=50,hl=100,xlim=c(45,105),ylim=c(0,105))
    tradeoffplot3(B50,VY,"Prob. biomass above half BMSY (%)", "Prob. AAVY less than 15%",MSEobj@MPs[1:MSEobj@nMPs],vl=80,hl=50,xlim=c(75,105),ylim=c(45,105))
  }
  
  #if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  #if(!is.na(nam) & !is.character(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
  #if(!is.na(nam) & is.character(nam))mtext(nam,3,outer=T,line=0.3,font=2)
  
  
  temp<-data.frame(PNOF,B50,LTY,VY)
  row.names(temp)<-MSEobj@MPs[1:MSEobj@nMPs]
  temp
  
}

# Projection plot - plots F/FMSY and B/BMSY for projection years and all MPs in 
# MSEobj.  
Pplot<-function(MSEobj,nam=NA){
  
  FMSYr<-quantile(MSEobj@F_FMSY,c(0.001,0.90),na.rm=T)
  BMSYr<-quantile(MSEobj@B_BMSY,c(0.001,0.975),na.rm=T)
  
  colsse<-rainbow(100,start=0,end=0.36)[1:100]
  colB<-rep(colsse[100],ceiling(BMSYr[2]*100))
  colB[1:100]<-colsse
  colB<-makeTransparent(colB,60)
  colsse<-rainbow(200,start=0,end=0.36)[200:1]
  colF<-rep(colsse[200],ceiling(FMSYr[2]*100))
  colF[1:200]<-colsse
  colF<-makeTransparent(colF,60)
  
  Yd<-rep(NA,MSEobj@nMPs)
  P10<-rep(NA,MSEobj@nMPs)
  P50<-rep(NA,MSEobj@nMPs)
  P100<-rep(NA,MSEobj@nMPs)
  POF<-rep(NA,MSEobj@nMPs)
  yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  RefYd<-MSEobj@OM$RefY
  
  for(mm in 1:MSEobj@nMPs){
    Yd[mm]<-round(mean(apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd,na.rm=T)*100,1)
    #cbind(MSEobj@C[,mm,yind],unlist(MSEobj@OM$MSY))
    POF[mm]<-round(sum(MSEobj@F_FMSY[,mm,]>1,na.rm=T)/prod(dim(MSEobj@F_FMSY[,mm,]),na.rm=T)*100,1)
    P10[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    P50[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<0.5,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    P100[mm]<-round(sum(MSEobj@B_BMSY[,mm,]<1,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
  }
    
  nr<-ceiling(MSEobj@nMPs/8)
  nc<-ceiling(MSEobj@nMPs/nr)
  nr<-nr*2
  MSEcols<-c('red','green','blue','orange','brown','purple','dark grey','violet','dark red','pink','dark blue','grey')
  temp<-array(0,c(nr*2+(nr/2-1),nc*2))
  i<-0
  for(c in 1:nc){
    for(r in 1:nr){
      i<-i+1
      temp[(ceiling(r/2)-1)+(1:2)+(r-1)*2,(1:2)+(c-1)*2]<-((c-1)*nr)+r
    }
  }
  par(mfcol=c(nr,nc),mai=c(0.2,0.35,0.3,0.01),omi=c(0.5,0.4,0.4,0.05))
  layout(temp)
  #dev.new2(width=nc*3,height=nr*3)
  #
  lwdy<-2.5

  for(mm in 1:MSEobj@nMPs){
    plot(MSEobj@F_FMSY[1,mm,],ylim=FMSYr,col=colF[ceiling(mean(MSEobj@F_FMSY[1,mm,],na.rm=T)*100)],type='l',lwd=lwdy)
    for(i in 1:MSEobj@nsim)lines(MSEobj@F_FMSY[i,mm,],col=colF[ceiling(mean(MSEobj@F_FMSY[i,mm,],na.rm=T)*100)],lwd=lwdy)
    abline(h=100,col="grey",lwd=3)
    mtext(MSEobj@MPs[mm],3,outer=F,line=0.6)
    legend('topright',c(paste(POF[mm],"% POF",sep=""),
                 paste(Yd[mm],"% FMSY yield",sep="")),bty='n',cex=0.8)
    if(mm%in%(1:(nr/2)))mtext("F/FMSY",2,line=2.5,outer=F)
    abline(h=1,col=makeTransparent("grey",30),lwd=2.5)
  
    plot(MSEobj@B_BMSY[1,mm,],ylim=BMSYr,col=colB[ceiling(MSEobj@B_BMSY[1,mm,MSEobj@proyears]*100)],type='l',lwd=lwdy)
    for(i in 1:MSEobj@nsim)lines(MSEobj@B_BMSY[i,mm,],col=colB[ceiling(MSEobj@B_BMSY[i,mm,MSEobj@proyears]*100)],lwd=lwdy)
    abline(h=100,col="grey",lwd=3)
    legend('topright',c(paste(P100[mm],"% < BMSY",sep=""),
                 paste(P50[mm],"% < 0.5BMSY",sep=""),
                 paste(P10[mm],"% < 0.1BMSY",sep="")),bty='n',cex=0.8)
    if(mm%in%(1:(nr/2)))mtext("B/BMSY",2,line=2.5,outer=F)
    abline(h=1,col=makeTransparent("grey",30),lwd=2.5)
  
  }
  mtext("Projection year",1,outer=T,line=1.2)
  if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.3,font=2)
  if(!is.na(nam))mtext(MSEobj@Name,3,outer=T,line=0.3,font=2)
}

# Kobe plot
Kplot<-function(MSEobj,maxsim=60,nam=NA){
  nr<-floor((MSEobj@nMPs)^0.5)
  nc<-ceiling((MSEobj@nMPs)/nr)
    
  FMSYr<-quantile(MSEobj@F_FMSY,c(0.001,0.90),na.rm=T)
  BMSYr<-quantile(MSEobj@B_BMSY,c(0.001,0.975),na.rm=T)
    
  #dev.new2(width=nc*3,height=nr*3.6)
  par(mfrow=c(nr,nc),mai=c(0.45,0.45,0.45,0.01),omi=c(0.45,0.3,0.35,0.01))
  
  colsse<-rainbow(MSEobj@proyears,start=0.63,end=0.95)[1:MSEobj@proyears]
  colsse<-makeTransparent(colsse,95)
  
  for(mm in 1:MSEobj@nMPs){
    plot(c(MSEobj@B_BMSY[1,mm,1],MSEobj@B_BMSY[1,mm,2]),
         c(MSEobj@F_FMSY[1,mm,1],MSEobj@F_FMSY[1,mm,2]),xlim=BMSYr,ylim=FMSYr,
         col=colsse[1],type='l')
    
    OO<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]<1&MSEobj@F_FMSY[,mm,MSEobj@proyears]>1,na.rm=T)/MSEobj@nsim*100,1)
    OU<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]>1&MSEobj@F_FMSY[,mm,MSEobj@proyears]>1,na.rm=T)/MSEobj@nsim*100,1)
    UO<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]<1&MSEobj@F_FMSY[,mm,MSEobj@proyears]<1,na.rm=T)/MSEobj@nsim*100,1)
    UU<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]>1&MSEobj@F_FMSY[,mm,MSEobj@proyears]<1,na.rm=T)/MSEobj@nsim*100,1)
    
    #alp<-80
    #polygon(c(1,-1000,-1000,1),c(1,1,1000,1000),col=makeTransparent("orange",alp),border=makeTransparent("orange",alp))
    #polygon(c(1,1000,1000,1),c(1,1,1000,1000),col=makeTransparent("yellow",alp),border=makeTransparent("yellow",alp))
    #polygon(c(1,-1000,-1000,1),c(1,1,-1000,-1000),col=makeTransparent("yellow",alp),border=makeTransparent("yellow",alp))
    #polygon(c(1,1000,1000,1),c(1,1,-1000,-1000),col=makeTransparent("green",alp),border=makeTransparent("yellow",alp))
    
    
    abline(h=1,col="grey",lwd=3)
    abline(v=1,col="grey",lwd=3)
    #abline(v=c(0.1,0.5),col="grey",lwd=2)
    rng<-1:min(maxsim,MSEobj@nsim)
    for(i in rng){
      for(y in 1:(MSEobj@proyears-1)){
        lines(c(MSEobj@B_BMSY[i,mm,y],MSEobj@B_BMSY[i,mm,y+1]),
              c(MSEobj@F_FMSY[i,mm,y],MSEobj@F_FMSY[i,mm,y+1]),
              col=colsse[y],lwd=1.6)
      }
    }
    
    points(MSEobj@B_BMSY[rng,mm,1],MSEobj@F_FMSY[rng,mm,1],pch=19,cex=0.8,col=colsse[1])
    points(MSEobj@B_BMSY[rng,mm,MSEobj@proyears],MSEobj@F_FMSY[rng,mm,MSEobj@proyears],pch=19,cex=0.8,col=colsse[MSEobj@proyears])
    
    if(mm==1)legend('right',c("Start","End"),bty='n',text.col=c(colsse[1],colsse[MSEobj@proyears]),pch=19,col=c(colsse[1],colsse[MSEobj@proyears]))
    legend('topleft',paste(OO,"%",sep=""),bty='n',text.font=2)
    legend('topright',paste(OU,"%",sep=""),bty='n',text.font=2)
    legend('bottomleft',paste(UO,"%",sep=""),bty='n',text.font=2)
    legend('bottomright',paste(UU,"%",sep=""),bty='n',text.font=2)
    
    mtext(MSEobj@MPs[mm],3,line=0.45)
  }
  mtext("B/BMSY",1,outer=T,line=1.4)
  mtext("F/FMSY",2,outer=T,line=0.2)
  if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.25,font=2)
  if(!is.na(nam))mtext(MSEobj@Name,3,outer=T,line=0.25,font=2)
}

# A simulation by simulation approach to plotting results
comp<-function(MSEobj,MPs=NA){
  
  if(is.na(MPs))MPs<-MSEobj@MPs
  notm<-MPs[!(MPs%in%MSEobj@MPs)]
  canm<-MPs[MPs%in%MSEobj@MPs]
  if(length(notm)>0)print(paste("Methods",paste(notm,collapse=", "),"were not carried out in MSE",deparse(substitute(MSEobj)),sep=" "))
  
  if(length(canm)==0)stop(paste('None of the methods you specified were carried out in the MSE', deparse(substitute(MSEobj)),sep=""))
  
  if(length(canm)>4){
    print(paste('A maximum of four methods can be compared at once. Plotting first four:',paste(canm[1:4],collapse=", "),sep=" "))
    canm<-canm[1:4]
  } 
  
  mind<-match(canm,MSEobj@MPs)
  nm<-length(mind)
  nsim<-MSEobj@nsim
  proyears<-MSEobj@proyears
  
  Yd<-array(NA,c(nm,nsim))
  P10<-array(NA,c(nm,nsim))
  P50<-array(NA,c(nm,nsim))
  P100<-array(NA,c(nm,nsim))
  POF<-array(NA,c(nm,nsim))
  yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  RefYd<-MSEobj@OM$RefY
  
  for(m in 1:nm){
    mm<-mind[m]
    Yd[m,]<-round(apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd*100,1)
    POF[m,]<-round(apply(MSEobj@F_FMSY[,mm,]>1,1,sum,na.rm=T)/proyears*100,1)
    P10[m,]<-round(apply(MSEobj@B_BMSY[,mm,]<0.1,1,sum,na.rm=T)/proyears*100,1)
    P50[m,]<-round(apply(MSEobj@B_BMSY[,mm,]<0.5,1,sum,na.rm=T)/proyears*100,1)
    P100[m,]<-round(apply(MSEobj@B_BMSY[,mm,]<1,1,sum,na.rm=T)/proyears*100,1)
  }
  
  MSEcols<-c('red','green','blue','orange')
  
  #dev.new2(width=7,height=7)
  par(mfrow=c(2,2),mai=c(0.85,0.7,0.1,0.1),omi=rep(0.01,4))
  
  tradeoffplot2(POF,Yd,"Prob. of overfishing (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=NA)
  tradeoffplot2(P100,Yd,"Prob. biomass < BMSY (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=canm)
  tradeoffplot2(P50,Yd,"Prob. biomass < 0.5BMSY (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=NA)
  tradeoffplot2(P10,Yd,"Prob. biomass < 0.1BMSY (%)", "Relative yield",vl=50,hl=100,coly=MSEcols,leg=NA)
  
}


# Supporting Functions for Plots
tradeoffplot<-function(x,y,xlab,ylab,labs,cex,vl,hl){
   adjj<-c(0.7,1.3)
   XLim <- c(min(c(-10, min(x,na.rm=T)*adjj)), max(c(max(x,na.rm=T)*adjj, 110)))
   YLim <- c(min(c(-10, min(y,na.rm=T)*adjj)), max(c(max(y,na.rm=T)*adjj, 110)))
   coly<-rep(c('#0000ff95','#ff000095','#20ff1095'),50)[1:length(labs)]
   coly[labs%in%c("AvC","curE","FMSYref")]<-'black'
   # plot(NA,xlim=range(x,na.rm=T)*adjj,ylim=range(y,na.rm=T)*adjj,xlab=xlab,ylab=ylab)
   plot(NA,xlim=XLim,ylim=YLim,xlab=xlab,ylab=ylab)
   abline(v=vl,col="#99999940",lwd=2)
   abline(h=hl,col="#99999940",lwd=2)
   text(x,y,labs,font=2,col=coly,cex=0.9)
}

tradeoffplot2<-function(x,y,xlab,ylab,cex=1,vl,hl,coly,leg){
  adjj<-c(0.7,1.3)
  plot(NA,xlim=range(x,na.rm=T)*adjj,ylim=range(y,na.rm=T)*adjj,xlab=xlab,ylab=ylab)
  abline(v=vl,col="grey",lwd=2)
  abline(h=hl,col="grey",lwd=2)
  for(m in 1:nrow(x))points(x[m,],y[m,],col=makeTransparent(coly[m],50),pch=19,cex=cex)
  if(!is.na(leg[1]))legend('topright',legend=leg,text.col=coly,bty='n')
}

tradeoffplot3<-function(x,y,xlab,ylab,labs,cex,vl,hl,xlim,ylim){
  coly<-rep(c('#0000ff95','#ff000095','#20ff1095'),10)
  plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
  abline(v=vl,col="#99999940",lwd=2)
  abline(h=hl,col="#99999940",lwd=2)
  text(x,y,labs,font=2,col=coly,cex=1)
}

makeTransparent<-function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


# Value of information analysis ------------------------------------------------
# Value of information
VOI<-function(MSEobj,ncomp=6,nbins=8,maxrow=8,Ut=NA,Utnam="Utility"){
  objnam<-deparse(substitute(MSEobj))
  nsim<-MSEobj@nsim
 
  if(is.na(Ut[1])){
    Ut<-array(NA,c(nsim,MSEobj@nMPs))
    yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
    RefYd<-MSEobj@OM$RefY
  
    for(mm in 1:MSEobj@nMPs){
      Ut[,mm]<-apply(MSEobj@C[,mm,yind],1,mean,na.rm=T)/RefYd*100
    #POF[,mm]<-apply(MSEobj@F_FMSY[,mm,]>1,1,sum)/MSEobj@proyears
    #P10[,mm]<-apply(MSEobj@B_BMSY[,mm,]<0.1,1,sum)/MSEobj@proyears
    }
    Utnam<-"Long-term yield relative to MSY (%)"
  }

  MPs<-MSEobj@MPs
  nMPs<-MSEobj@nMPs
 
  onlycor<-c("RefY","A","MSY","Linf","t0","OFLreal","Spat_targ")
  vargood<-(apply(MSEobj@OM,2,sd)/(apply(MSEobj@OM,2,mean)^2)^0.5)>0.005
  # MSEobj@OM<- MSEobj@OM[,(!names(MSEobj@OM)%in%onlycor)&vargood]
  MSEobj@OM[,which((!names(MSEobj@OM)%in%onlycor)&vargood)]  
  OMp<-apply(MSEobj@OM,2,quantile,p=seq(0,1,length.out=nbins+1))
  Obsp<-apply(MSEobj@Obs,2,quantile,p=seq(0,1,length.out=nbins+1))
  OMv<-array(NA,c(nMPs,ncol(MSEobj@OM),nbins))
  Obsv<-array(NA,c(nMPs,ncol(MSEobj@Obs),nbins))
  
  for(mm in 1:nMPs){
    for(j in 1:nbins){
      for(i in 1:ncol(MSEobj@OM)){
        cond<-MSEobj@OM[,i]>OMp[j,i]&MSEobj@OM[,i]<OMp[j+1,i]
        OMv[mm,i,j]<-mean(Ut[cond,mm],na.rm=T)
      }
      for(i in 1:ncol(MSEobj@Obs)){
        cond<-MSEobj@Obs[,i]>Obsp[j,i]&MSEobj@Obs[,i]<Obsp[j+1,i]
        Obsv[mm,i,j]<-mean(Ut[cond,mm],na.rm=T)
      }
    }
  }
  
  # -- Operating model variables
  OMs<-apply(OMv,1:2,sd,na.rm=T)
  OMstr<-array("",c(nMPs*2,ncomp+1))
   
  for(mm in 1:nMPs){
    ind<-order(OMs[mm,],decreasing=T)[1:ncomp]
    OMstr[1+(mm-1)*2,1]<-MPs[mm]
    OMstr[1+(mm-1)*2,2:(1+ncomp)]<-names(MSEobj@OM[ind])
    OMstr[2+(mm-1)*2,2:(1+ncomp)]<-round(OMs[mm,ind],2)
  }
  OMstr<-data.frame(OMstr)
  names(OMstr)<-c("MP",1:ncomp)
  
  
  # -- Observation model variables
  slots<-c( "Cat",  "Cat","AvC",  "AvC","CAA",      "CAA",    "CAL",      "CAL",    "Ind","Dep",  "Dep", "Dt",   "Dt", "Mort", "FMSY_M",    "BMSY_B0",     "L50",      "L95",    "LFC",    "LFS",    "Abun",  "Abun","vbK",  "vbt0",  "vbLinf",  "Steep","Iref",    "Cref",    "Bref")
  Obsnam<-c("Cbias","Csd","Cbias","Csd","CAA_nsamp","CAA_ESS","CAL_nsamp","CAL_ESS","Isd","Dbias","Derr","Dbias","Derr","Mbias","FMSY_Mbias","BMSY_B0bias", "lenMbias","lenMbias","LFCbias","LFSbias","Abias","Aerr","Kbias","t0bias","Linfbias","hbias","Irefbias","Crefbias","Brefbias")
  Obss<-apply(Obsv,1:2,sd,na.rm=T)
  Obsstr<-array("",c(nMPs*2,ncomp+1))
  for(mm in 1:nMPs){
    relobs<-Obsnam[slots%in%unlist(strsplit(Required(MPs[mm])[,2],split=", "))]
    ind<-(1:ncol(MSEobj@Obs))[match(relobs,names(MSEobj@Obs))]
    pos<-names(MSEobj@Obs)[ind]# possible observation processes
    maxy<-min(max(1,length(pos)),ncomp,na.rm=T)
    ind2<-order(Obss[mm,ind],decreasing=T)[1:maxy]
    Obsstr[1+(mm-1)*2,1]<-MPs[mm]
    Obsstr[1+(mm-1)*2,2:(1+maxy)]<-pos[ind2]
    Obsstr[2+(mm-1)*2,2:(1+maxy)]<-round(Obss[mm,ind][ind2],2)
  }
  Obsstr<-data.frame(Obsstr)
  names(Obsstr)<-c("MP",1:ncomp)

  ncols<-40
  #colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),95)[ncols:1]
  colsse<-makeTransparent(rainbow(ncols,start=0,end=0.36),90)[ncols:1]
  minsd<-0
  maxsd<-max(OMs,na.rm=T)
  coly<-ceiling(OMs/maxsd*ncols)
 
  # Operating model variables
  mbyp <- split(1:nMPs, ceiling(1:nMPs/maxrow))
  ylimy=c(0,max(OMv,na.rm=T)*1.2)              
  
  
  for(pp in 1:length(mbyp)){
    
    par(mfrow=c(length(mbyp[[pp]]),ncomp),mai=c(0.15,0.1,0.15,0.05),omi=c(0.1,0.9,0.3,0.05))
   
    for(mm in mbyp[[pp]]){
      for(cc in 1:ncomp){
        rind<-(mm-1)*2+1
        y<-Ut[,mm]
        cind<-match(OMstr[rind,1+cc],names(MSEobj@OM))
        x<-MSEobj@OM[,cind]
        plot(x,y,col="white",axes=F,ylim=ylimy)
        axis(1,pretty(OMp[,cind]),pretty(OMp[,cind]),cex.axis=0.8,padj=-1.5)
        abline(v=OMp[,cind],col="#99999960")
        points(x,y,col=colsse[coly[mm,cind]],pch=19,cex=0.8)
        x2<-(OMp[1:nbins,cind]+OMp[2:(nbins+1),cind])/2
        y2<-OMv[mm,cind,]
        lines(x2,y2)
        legend('bottomright',legend=round(OMs[mm,cind],2),bty='n',cex=0.8)
        legend('topleft',legend=OMstr[rind,1+cc],bty='n',cex=0.85)
        if(cc==1){ 
          mtext(MPs[mm],2,font=2,outer=F,cex=0.8,line=2)
          ytick<-pretty(seq(ylimy[1],ylimy[2]*1.3,length.out=10))
          axis(2,ytick,ytick,cex.axis=0.8)
        } # only first column
      } # parameters (columns)
    } # MPs (rows)
    
    mtext(Utnam,2,outer=T,cex=0.9,line=3.5)
    mtext(paste("Operating model parameters: ",objnam,"@OM",sep=""),3,outer=T,font=2,cex=0.9)
    
  } # Plots
  
  # Observation model values
  
  ylimy=c(0,max(Obsv,na.rm=T)*1.2)              
  minsd<-0
  maxsd<-max(Obss)
  coly<-ceiling(Obss/maxsd*ncols)
  
  if(sum(is.na(Obsstr)|Obsstr=="")<(ncomp+1)*nMPs*2-nMPs){ # only if there is data to plot
  
  for(pp in 1:length(mbyp)){
    
    par(mfrow=c(length(mbyp[[pp]]),ncomp),mai=c(0.15,0.1,0.15,0.05),omi=c(0.1,0.9,0.3,0.05))
    
    for(mm in mbyp[[pp]]){
      rind<-(mm-1)*2+1
      npres<-sum(Obsstr[rind+1,]!="")
      for(cc in 1:ncomp){
        if(!is.na(npres)&cc<(npres+1)){
          y<-Ut[,mm]
          cind<-match(Obsstr[rind,1+cc],names(MSEobj@Obs))
          x<-MSEobj@Obs[,cind]
          plot(x,y,col="white",axes=F,ylim=ylimy)
          axis(1,pretty(Obsp[,cind]),pretty(Obsp[,cind]),cex.axis=0.8,padj=-2)
          abline(v=Obsp[,cind],col="#99999960")
          points(x,y,col=colsse[coly[mm,cind]],pch=19,cex=0.8)
          x2<-(Obsp[1:nbins,cind]+Obsp[2:(nbins+1),cind])/2
          y2<-Obsv[mm,cind,]
          lines(x2,y2)
          legend('bottomright',legend=round(Obss[mm,cind],2),bty='n',cex=0.8)
          legend('topleft',legend=Obsstr[rind,1+cc],bty='n',cex=0.75)
          if(cc==1){ 
            mtext(MPs[mm],2,font=2,outer=F,cex=0.6,line=2)
            ytick<-pretty(seq(ylimy[1],ylimy[2]*1.3,length.out=10))
            axis(2,ytick,ytick,cex.axis=0.8)
          } # only first column
        }else{
          plot(0,type='n',axes=FALSE,ann=FALSE)
          if(cc==1){ 
            mtext(MPs[mm],2,font=2,outer=F,cex=0.6,line=2)
          } # only first column
        }  
      } # parameters (columns)
    } # MPs (rows)
    
    mtext(Utnam,2,outer=T,cex=0.9,line=3.5)
    mtext(paste("Observation model parameters: ",objnam,"@Obs",sep=""),3,outer=T,font=2,cex=0.9)
   
  } # Plots
  } # if there is data to plot
  
  list(OMstr,Obsstr)
  
} # VOI

# Manipulation of MSE Object ---------------------------------------------------
# Subset the MSE object by particular MPs (either MP number or name), 
#  or particular simulations
Sub <- function(MSEobj, MPs=NULL, sims=NULL) {
  Class <- class(MPs)
  if(Class == "NULL") subMPs <- MSEobj@MPs
  if(Class == "integer" | Class == "numeric") subMPs <- MSEobj@MPs[as.integer(MPs)]
  if(Class == "character") subMPs <- MPs
  SubMPs <- which(MSEobj@MPs %in% subMPs)
  not <- (subMPs %in% MSEobj@MPs) # Check for MPs misspelled
  ind <- which(not == FALSE)
  newMPs <- MSEobj@MPs[SubMPs]
  if (length(SubMPs) < 1) stop("MPs not found")
  if (length(ind > 0)) {
    message("Warning: MPs not found - ", paste0(subMPs[ind], " "))
	message("Subsetting by MPs: ", paste0(newMPs, " "))
  }
  
  ClassSims <- class(sims)
  if (ClassSims == "NULL") SubIts <- 1:MSEobj@nsim
  if (ClassSims == "integer" | ClassSims == "numeric") {
    # sims <- 1:min(MSEobj@nsim, max(sims))
	SubIts <- as.integer(sims)
  }	
  if (ClassSims == "logical") SubIts <- which(sims)

  SubF <- MSEobj@F_FMSY[SubIts,SubMPs,]
  SubB <- MSEobj@B_BMSY[SubIts,SubMPs,]
  OutOM <- MSEobj@OM[SubIts,]
  
  SubResults <- new('MSE',Name=MSEobj@Name, nyears=MSEobj@nyears, proyears=MSEobj@proyears, nMPs=length(SubMPs),
	MPs=newMPs, nsim=length(SubIts), OMtable=OutOM, Obs=MSEobj@Obs[SubIts,], B_BMSYa=SubB, F_FMSYa=SubF, Ba=MSEobj@B[SubIts,SubMPs,], 
	FMa=MSEobj@FM[SubIts,SubMPs,], Ca=MSEobj@C[SubIts,SubMPs,], TACa=MSEobj@TAC[SubIts,SubMPs,], SSB_hist=MSEobj@SSB_hist[SubIts,,,],
	CB_hist=MSEobj@CB_hist[SubIts,,,], FM_hist=MSEobj@FM_hist[SubIts,,,])
  
 return(SubResults)
}


# Evaluate Peformance of MPs ---------------------------------------------------
# Function examines how consistently an MP outperforms another. 
DOM <- function(MSEobj, MPtg=NA) {
  if (any(is.na(MPtg))) MPtg <- MSEobj@MPs 
  proyears<-MSEobj@proyears
  nMP <- MSEobj@nMPs
  ind <- which(MSEobj@MPs %in%  MPtg)
  MPr <- which(!(MSEobj@MPs %in%  MPtg))
  yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
  y1<-1:(MSEobj@proyears-1)
  y2<-2:MSEobj@proyears
  Mat <- matrix(0, nrow=length(MPtg), ncol=nMP)
  rownames(Mat) <- MPtg
  colnames(Mat) <- MSEobj@MPs
  POF <- P100 <- YieldMat <- IAVmat <- Mat
  for (X in 1:length(MPtg)) {
    # Overfishing (F > FMSY)
    ind1 <- as.matrix(expand.grid(1:nsim, ind[X], 1:proyears))
    ind2 <- as.matrix(expand.grid(1:nsim,  1:nMP, 1:proyears))
    t1 <- apply(array(MSEobj@F_FMSY[ind1] > 1, dim=c(nsim, 1, proyears)), c(1,2), sum, na.rm=TRUE) 
    t2 <- apply(array(MSEobj@F_FMSY[ind2] > 1, dim=c(nsim, nMP, proyears)), c(1,2), sum, na.rm=TRUE)
    POF[X,] <- round(apply(matrix(rep(t1, nMP), nrow=nsim) < t2, 2, sum) / nsim * 100, 0)
	# B < BMSY
    t1 <- apply(array(MSEobj@B_BMSY[ind1] < 1, dim=c(nsim, 1, proyears)), c(1,2), sum, na.rm=TRUE) 
    t2 <- apply(array(MSEobj@B_BMSY[ind2] < 1, dim=c(nsim, nMP, proyears)), c(1,2), sum, na.rm=TRUE) 
    P100[X,] <- round(apply(matrix(rep(t1, nMP), nrow=nsim) < t2, 2, sum, na.rm=TRUE) / nsim * 100, 0)
    # Relative yield in last 5 years
    ind1 <- as.matrix(expand.grid(1:nsim, ind[X], yind))
    ind2 <- as.matrix(expand.grid(1:nsim,  1:nMP, yind))
    t1 <- apply(array(MSEobj@C[ind1], dim=c(nsim, 1, length(yind))), c(1,2), sum, na.rm=TRUE)
    t2 <- apply(array(MSEobj@C[ind2], dim=c(nsim, nMP, length(yind))), c(1,2), sum, na.rm=TRUE)
    YieldMat[X,] <- round(apply(matrix(rep(t1, nMP), nrow=nsim) > t2, 2, sum, na.rm=TRUE) / nsim * 100, 0)
    # interannual variation in catch 
    ind1 <- as.matrix(expand.grid(1:nsim, ind[X], y1))
    ind2 <- as.matrix(expand.grid(1:nsim, ind[X], y2)) 
    AAVY1 <- apply(array(((MSEobj@C[ind1]-MSEobj@C[ind2])^2)^0.5, dim=c(nsim, 1, length(y1))),1,mean,na.rm=T)/apply(array(MSEobj@C[ind2], dim=c(nsim, 1, length(y1))),1,mean,na.rm=T) 
    ind1 <- as.matrix(expand.grid(1:nsim, 1:nMP, y1))
    ind2 <- as.matrix(expand.grid(1:nsim, 1:nMP, y2)) 
    AAVY2 <- apply(array(((MSEobj@C[ind1]-MSEobj@C[ind2])^2)^0.5, dim=c(nsim, nMP, length(y1))),c(1,2),mean,na.rm=T)/apply(array(MSEobj@C[ind2], dim=c(nsim, nMP, length(y1))),c(1,2),mean,na.rm=T)
    IAVmat[X,] <- round(apply(matrix(rep(AAVY1, nMP), nrow=nsim) < AAVY2, 2, sum, na.rm=TRUE) / nsim * 100, 0)
  }
  out <- list() 
  out$POF <- POF
  out$P100 <- P100
  out$Yd <- YieldMat
  out$AAVY <- IAVmat
  return(out)
}


# Trade-Off Plot Function ------------------------------------------------------
TradePlot <- function(MSEobj, XAxis=c("Overfishing", "Biomass:BMSY"), 
	YAxis=c("Long-term Yield", "AnnualVar"), XThresh=c(30, 80), YThresh=c(0,50),
	maxVar=15, BmsyRef=0.5, B0Ref=0.2, AvailMPs=NULL, ShowLabs=FALSE, 
	ShowCols=TRUE) {
  PMs <- c("Long-term Yield", "Short-term Yield", "Overfishing", "Biomass:BMSY",
	"Biomass:B0", "AnnualVar")
  # Error Checks 	
  if (prod(XAxis %in% PMs)!=1) {
    message("Available Performance Metrics")
    print(PMs)
    stop("Invalid XAxis Performance Metrics")
  }	
  if (prod(YAxis %in% PMs)!=1) {
    message("Available Performance Metrics")
    print(PMs)
    stop("Invalid YAxis Performance Metrics")
  }	
  if (length(XAxis) > 4) stop("Too many Performance Metrics (max 4)")
  if (length(YAxis) > 4) stop("Too many Performance Metrics (max 4)")
  if (length(XAxis) != length(YAxis)) stop("XAxis must be of same length as YAxis")
  if (length(XThresh) != length(XAxis) | length(YThresh) != length(XAxis)) 
	warning("Risk Threshold not same length as number of PMs")
  
  Yd<-rep(NA,MSEobj@nMPs)
  BMSYref<-rep(NA,MSEobj@nMPs)
  B0ref<-rep(NA,MSEobj@nMPs)
  PNOF<-rep(NA,MSEobj@nMPs)
  LTY<-rep(NA,MSEobj@nMPs)
  STY<-rep(NA,MSEobj@nMPs)
  VY<-rep(NA,MSEobj@nMPs)

  y1 <- 1:(MSEobj@proyears-1)
  y2 <- 2:MSEobj@proyears
  
  ystart<-1:5
  yend<-max(MSEobj@proyears-4,1):MSEobj@proyears
  
  RefYd<-MSEobj@OM$RefY
  if (maxVar < 1) maxVar <- maxVar * 100
  
  for(mm in 1:MSEobj@nMPs){  
    PNOF[mm]<-round(sum(MSEobj@F_FMSY[,mm,]<1,na.rm=T)/prod(dim(MSEobj@F_FMSY[,mm,]),na.rm=T)*100,1)
    BMSYref[mm]<-round(sum(MSEobj@B_BMSY[,mm,]>BmsyRef,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
	B0ref[mm]<-round(sum((MSEobj@B[,mm,]/MSEobj@B[,mm,1] * MSEobj@OM$Depletion) >B0Ref,na.rm=T)/prod(dim(MSEobj@B_BMSY[,mm,]))*100,1)
    # LTY[mm]<-round(sum(MSEobj@C[,mm,yend]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(yend)),3)*100
	# STY[mm]<-round(sum(MSEobj@C[,mm,ystart]/RefYd>0.5,na.rm=T)/(MSEobj@nsim*length(ystart)),3)*100
	LTY[mm]<-round(mean(apply(MSEobj@C[,mm,yend],1,mean,na.rm=T)/RefYd,na.rm=T)*100,1)
	STY[mm]<-round(mean(apply(MSEobj@C[,mm,ystart],1,mean,na.rm=T)/RefYd,na.rm=T)*100,1)
    AAVY<-apply((((MSEobj@C[,mm,y1]-MSEobj@C[,mm,y2])/MSEobj@C[,mm,y2])^2)^0.5,1,mean,na.rm=T) 
    VY[mm]<-round(sum(AAVY<(maxVar/100),na.rm=T)/MSEobj@nsim,3)*100
  }
  
  for (xx in seq_along(XAxis)) {
    name <- paste0("X", xx)
	name1 <- paste0("XLab", xx)
    assign(name, GetStat(XAxis[xx], LTY, STY, PNOF, BMSYref, B0ref, VY))
	assign(name1, StatLab(XAxis[xx], maxVar, BmsyRef, B0Ref))
	name <- paste0("Y", xx)
	name1 <- paste0("YLab", xx)
	assign(name, GetStat(YAxis[xx], LTY, STY, PNOF, BMSYref, B0ref, VY))
	assign(name1, StatLab(YAxis[xx], maxVar, BmsyRef, B0Ref))
  }
  
  Nplot <- length(XAxis)
  if (Nplot == 1) par(mfrow=c(1,1), mar=c(4,4.5,1,1), oma=c(1,1,0,0))
  if (Nplot == 2) par(mfrow=c(1,2), mar=c(4,4.5,1,1), oma=c(1,1,0,0))
  if (Nplot == 3) par(mfrow=c(1,3), mar=c(4,4.5,1,1), oma=c(1,1,0,0))
  if (Nplot == 4) par(mfrow=c(2,2), mar=c(4,4.5,1,1), oma=c(1,1,0,0))
  
  OutList <- list()
  for (xx in seq_along(XAxis)) {
    Xname <- paste0("X", xx)
	XLab <- paste0("XLab", xx)
	Yname <- paste0("Y", xx)
	YLab <- paste0("YLab", xx)
    rr <- tradeoffplot4(x=get(Xname), y=get(Yname), get(XLab), get(YLab), 
		labs=MSEobj@MPs[1:MSEobj@nMPs],vl=XThresh[xx],hl=YThresh[xx], 
		ShowLabs=ShowLabs,  ShowCols=ShowCols, AvailMPs=AvailMPs)
	
	labs <- MSEobj@MPs[1:MSEobj@nMPs]
	ind <- which(labs %in% rr)
    tempDF <- data.frame(MP=rr, X=get(Xname)[ind], Y=get(Yname)[ind])
	Dist <- NULL # calculate distance from corner
    for (X in 1:length(tempDF[,2])) Dist[X] <- euc.dist(c(tempDF[X,2], tempDF[X,3]), c(100, 100))
	tempDF <- tempDF[order(Dist),]
	rownames(tempDF) <- 1:nrow(tempDF)
	OutList[[xx]] <- tempDF
  }
 
  print(OutList)
  invisible(OutList)
  
}

# Supporting functions 
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

GetStat <- function(PM, LTY, STY, PNOF, BMSYref, B0ref, VY) {
  switch(PM,
    "Long-term Yield" = LTY,
	"Short-term Yield" = STY,
	"Overfishing" = PNOF,
	"Biomass:BMSY" = BMSYref,
	"Biomass:B0" = B0ref,
    "AnnualVar" = VY)
}

StatLab <- function(PM, maxVar, BmsyRef, B0Ref) {
  switch(PM,
    "Long-term Yield" = "Long-term Yield",
	"Short-term Yield" = "Short-term Yield",
	"Overfishing" = "Prob. of Not Overfishing (%)",
	"Biomass:BMSY" = paste0("Prob. Biomass >",BmsyRef, "BMSY (%)"),
	"Biomass:B0" = paste0("Prob. Biomass >",B0Ref, "B0 (%)"),,
    "AnnualVar" = paste0("Prob. AAVY <", maxVar, "%")
	)
}

tradeoffplot4<-function(x,y,xlab,ylab,labs,cex,vl,hl, 
	ShowLabs=FALSE,  ShowCols=FALSE, AvailMPs=NULL){
   adjj<-c(0.9,1.1)
   XLim <- c(min(c(-10, min(x,na.rm=T)*adjj)), max(c(max(x,na.rm=T)*adjj, 110)))
   YLim <- c(min(c(-10, min(y,na.rm=T)*adjj)), max(c(max(y,na.rm=T)*adjj, 110)))
   
   # Which MPs meet minimum PMs 
   ind <- which(x >= vl & y >=hl)
   coly <- rep("darkgray", length(labs)) 
   coly[ind] <- "black" 
   # coly[labs%in%c("AvC","curE","FMSYref")]<-'black'
   
   Pch <- rep(21, length(labs))
   # Pch[labs%in%c("AvC","curE","FMSYref")] <- 17
   coly[grep("FMSY", labs)]<-'black'
   Pch[grep("FMSY", labs)] <- 24
   if (!is.null(AvailMPs)) Pch[labs%in%AvailMPs] <- 21
   if (!is.null(AvailMPs)) coly[labs%in%AvailMPs & (x >= vl & y >=hl)] <- "green"
   # coly<-rep(c('#0000ff95','#ff000095','#20ff1095'),50)[1:length(labs)]

   plot(NA,xlim=XLim,ylim=YLim,xlab=xlab,ylab=ylab, bty="l", las=1)
   abline(v=vl,col="#99999940",lwd=2)
   abline(h=hl,col="#99999940",lwd=2)
   
   Alpha <- 30
   # polygons 
   LeftCol <- rgb(red=255, green=0, blue=0, alpha=Alpha, names = NULL, 
	maxColorValue = 255)
   RightCol <- rgb(red=0, green=255, blue=0, alpha=Alpha, names = NULL, 
	maxColorValue = 255)   

   if(ShowCols) {
     polygon(x=c(0, vl,  vl, 0), y=c(0, 0, hl, hl), col=LeftCol, border=NA)
     polygon(x=c(0, vl,  vl, 0), y=c(0, 0, 100, 100), col=LeftCol, border=NA)
     polygon(x=c(vl,  100, 100, vl), y=c(0, 0, 100, 100), col=RightCol, border=NA)
     polygon(x=c(vl, 100,  100, vl), y=c(hl, hl, 100, 100), col=RightCol, border=NA)
    }
   
    Cex <- 1.5
   if(!ShowLabs) points(x,y, bg=coly, pch=Pch, cex=Cex, col="black" )
   if(ShowLabs) text(x,y,labs,font=2,col=coly,cex=1)
   # if(IdPoints) {
    # message("Click points on plot to display MP name")
	# message("Click Stop to finish")
	# flush.console()
	# identify(x,y, labels=labs)
   # }	
   
   labs[ind]
   
}


