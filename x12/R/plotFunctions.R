plot.x12work<-function(x,plots=c(1:9),...){
#plots 1: Original
#plots 2: Original Trend Adjusted
#plots 3: Log Original
#plots 4: Seasonal Factors
#plots 5: Seasonal Factors with SI Ratios
#plots 6: Spectrum Adjusted Orig
#plots 7: Spectrum Seasonal Adjusted
#plots 8: Spectrum Irregular
#plots 9: Spectrum Residulas
  par(ask=TRUE)
  if(x$seats)
    plots <- plots [apply(cbind(!plots==3,!plots==4,!plots==5),1,all)]
  if(any(plots==1)){
    plot_original(x)    
  }
  if(any(plots==2)){
    plot_original_seasonal_trend(x)
  }
  if(any(plots==3)){
    plot_original(x,log_transform=TRUE)
  }
  if(any(plots==4)){
    plot_seasonal_factors(x,SI_Ratios=FALSE)
  }
  if(any(plots==5)){
    plot_seasonal_factors(x)
  }
  if(any(plots==6)){
    plot_spectrum(x,which="original")
  }
  if(any(plots==7)){
    plot_spectrum(x,which="seasonaladj")
  }
  if(any(plots==8)){
    plot_spectrum(x,which="irregular")
  }
  if(any(plots==9)){
    plot_spectrum(x,which="residuals")
  }
  par(ask=FALSE)
}

plot_seasonal_factors <- function(out,SI_Ratios=TRUE,ylab="Value",xlab="",lwd_seasonal=1,
    col_seasonal="black",lwd_mean=1,col_mean="blue",col_siratio="darkgreen",col_replaced="red",
    cex_siratio=.9,cex_replaced=.9,SI_Ratios_replaced=TRUE,
    plot_legend=TRUE,legend_horiz=FALSE,legend_bty="o",
    
    ...){
  if(!SI_Ratios)
    v <- as.vector(out[["d10"]]) # Seasonal Factors
  else
    v <- as.vector(out[["d10"]])[1:length(out[["d8"]])] # Seasonal Factors without forecast
  f <- frequency(out[["d10"]])
  dif <- length(v)%%f
  if(dif>0)
    v[length(v)+(1:(f-dif))]<-NA
  out_matrix <- matrix(v,ncol=f,byrow=TRUE)
  if(f==12){
    lab <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  }else if(f==4){
    lab <- c("Qtr1","Qtr2","Qtr3","Qtr4")
  }else if(f==2){
    lab <- c("1st Half","2nd Half")
  }else{
    lab <- 1:f
  }
  if(SI_Ratios){
    main="Seasonal Factors by period and SI Ratios"
  }else{
    main="Seasonal Factors by period"
  }
  ylim <- c(min(v,na.rm=TRUE)*.95,max(v,na.rm=TRUE)*1.09)
  xlim <- c(0,f)
  plot(1,type="n",main=main,xlim=xlim,ylim=ylim,xaxt="n",ylab=ylab,xlab=xlab,cex=cex_siratio,...)
  axis(1,at=(1:f)-1/2,labels=lab)
  for(i in 0:(f)){    
    abline(v=i,col="grey")
  }
  if(SI_Ratios){
    vv <- as.vector(out[["d8"]]) #final unmodified SI Ratios
    dif <- length(vv)%%f
    if(dif>0)
      vv[length(vv)+(1:(f-dif))]<-NA
    out_matrix2 <- matrix(vv,ncol=f,byrow=TRUE)
    vvv <- as.vector(out[["d9"]]) # final replacement for SI Ratios
    dif <- length(vvv)%%f
    if(dif>0)
      vvv[length(vvv)+(1:(f-dif))]<-NA
    out_matrix3 <- matrix(vvv,ncol=f,byrow=TRUE)
  }
  for(i in 0:(f-1)){
    s <- seq(.1+i,(i+1)-.1,l=nrow(out_matrix))
    m <- mean(out_matrix[,i+1],na.rm=TRUE)
    points(rep(m,2)~c(s[1],s[length(s)]),type="l",col=col_mean,lwd=lwd_mean)
    points(out_matrix[,i+1]~s,type="l",col=col_seasonal,lwd=lwd_seasonal)
    if(SI_Ratios){
      points(out_matrix2[,i+1]~s,pch=20,cex=cex_siratio,col=col_siratio)
      if(SI_Ratios_replaced)
        points(out_matrix3[,i+1]~s,pch=20,cex=cex_replaced,col=col_replaced)
    }
  }
  if(plot_legend){
    if(SI_Ratios){
      if(SI_Ratios_replaced)
        legend(x=(f/2)-1,y=ylim[2],legend=c("Seasonal Factors","Mean","SI Ratio","Replaced SI Ratio"),
            col=c(col_seasonal,col_mean,col_siratio,col_replaced),pch=c(NA,NA,20,20),
            lty=c(1,1,NA,NA),bg="white",horiz=legend_horiz,bty=legend_bty)
      else
        legend(x=(f/2)-1,y=ylim[2],legend=c("Seasonal Factors","Mean","SI Ratio"),
            col=c(col_seasonal,col_mean,col_siratio),pch=c(NA,NA,20),
            lty=c(1,1,NA),bg="white",horiz=legend_horiz,bty=legend_bty)      
    }else
      legend(x=(f/2)-1,y=ylim[2],legend=c("Seasonal Factors","Mean"),col=c(col_seasonal,col_mean),
          lty=c(1,1),bg="white",horiz=legend_horiz,bty=legend_bty)
  }
}

plot_original <- function(out,ylab="Value",xlab="Date",
    main=if(!log_transform){"Original Series"}else{"Logs of the Original Series"},
    col="black",ytop=1,log_transform=FALSE,...){
  if(!log_transform)
    ts <- out[["a1"]]
  else
    ts <- log(out[["a1"]])
  plot(ts,ylim=c(min(ts,na.rm=TRUE),max(ts,na.rm=TRUE)*ytop),xlab=xlab,ylab=ylab,main=main,col=col,...)
}

plot_original_seasonal_trend <- function(out,ylab="Value",xlab="Date",
    main="Original Series, Seasonally Adjusted Series and Trend",
    col_original="black",col_seasonaladj="blue",col_trend="green",
    lwd_original=1,lwd_seasonaladj=1,lwd_trend=1,
    seasonaladj=TRUE,trend=TRUE,original=TRUE,plot_legend=TRUE,
    legend_horiz=TRUE,legend_bty="o",
    log_transform=FALSE,...){
  if(original)
    plot_original(out,ytop=1.1,col=col_original,main=main,xlab=xlab,ylab=ylab,lwd=lwd_original,
        log_transform=log_transform,...)
  else
    plot_original(out,ytop=1.1,col=col_original,main=main,xlab=xlab,ylab=ylab,
        lwd=lwd_original,log_transform=log_transform,type="n",...)
  text_leg <- vector()
  col_leg <- vector()
  if(original){
    text_leg <- "Original"
    col_leg <- c(col_original)
  }
  if(seasonaladj){
    ts_adj <- out[["d11"]]
    if(log_transform)
      ts_adj <- log(ts_adj)
    points(ts_adj,col=col_seasonaladj,type="l",lwd=lwd_seasonaladj)
    text_leg <- c(text_leg,"Seasonally Adjusted")
    col_leg <- c(col_leg,col_seasonaladj)
  }
  if(trend){
    ts_trend <- out[["d12"]]
    if(log_transform)
      ts_trend <- log(ts_trend)
    points(ts_trend,col=col_trend,type="l",lwd=lwd_trend)
    text_leg <- c(text_leg,"Trend")
    col_leg <- c(col_leg,col_trend)
  }
  lty <- rep(1,length(col_leg))
  if(plot_legend){
    if(!log_transform){
      legend(x=start(out[["a1"]])[1],y=max(out[["a1"]]*1.05,na.rm=TRUE)*1.05,lty=lty,legend=text_leg,
          col=col_leg,horiz=legend_horiz,bty=legend_bty)
    }else{
      legend(x=start(out[["a1"]])[1],y=log(max(out[["a1"]]*1.05,na.rm=TRUE))*1.05,lty=lty,
          legend=text_leg,col=col_leg,horiz=legend_horiz,bty=legend_bty)
    }
  }
}

plot_spectrum <- function(out,which="seasonaladj",xlab="Frequency",ylab="Decibels",
    main="default",
    col_bar="darkgrey",col_seasonal="red",col_td="blue",
    lwd_bar=4,lwd_seasonal=4,lwd_td=4,plot_legend=TRUE,
    legend_horiz=TRUE,legend_bty="o",
    ...){
  if(which=="seasonaladj")
    which <- "sp1"
  else if(which=="original")
    which <- "sp0"
  else if(which=="irregular")
    which <- "sp2"
  else if(which=="residuals")
    which <- "spr"
  if(main=="default"){
    if(which=="sp1")
      main <- "Spectrum of the Seasonally Adjusted Series"
    else if(which=="sp0")
      main <- "Spectrum of the Original Series"      
    else if(which=="sp2")
      main <- "Spectrum of the Irregular"
    else if(which=="spr")
      main <- "Spectrum of the RegARIMA Residuals"   
  }
  plot_spectrum_work(out[[which]]$frequency,out[[which]]$spectrum,xlab=xlab,ylab=ylab,
      main=main,
      col_bar=col_bar,col_seasonal=col_seasonal,col_td=col_td,
      lwd_bar=lwd_bar,lwd_seasonal=lwd_seasonal,lwd_td=lwd_td,plot_legend=plot_legend,
      legend_horiz=legend_horiz,legend_bty=legend_bty,
      ...)
}

plot_spectrum_work <- function(frequency,spectrum,xlab="Frequency",ylab="Decibels",
    f=12,main="default",highlight=TRUE,
    col_bar="darkgrey",col_seasonal="red",col_td="blue",
    lwd_bar=4,lwd_seasonal=4,lwd_td=4,plot_legend=TRUE,
    legend_horiz=TRUE,legend_bty="o",
    ...)
{
  gp<-par()
  for(i in c("cin","cra","csi","cxy","din")){
    gp <- gp[-which(names(gp)%in%i)]			
  }			
  tryCatch({
        par(mar = c(4, 4, 4, 2) + 0.1) 
        layout(matrix(c(rep(1,16),2,2),nrow=9,byrow=TRUE))#,heights = c(1, 6), respect = FALSE)
        plot(frequency,spectrum,type="n",xlab=xlab,ylab=ylab,main=main,col=col_bar,xaxt="n",...)
        
        #f <- 12#frequency(out[["a1"]])
#	abline(v=(1:(f/2))*1/f,col=col_seasonal,lwd=lwd_seasonal)
#	if(f==12)
#		abline(v=frequency[c(43,53)],col=col_td,lwd=lwd_td)
        coord <- par("usr")[3]
        
        if(highlight){	
          for(i in 1:length(frequency)){
            points(x=rep(frequency[i],2),y=c(spectrum[i],coord),type="l",col=col_bar,lwd=lwd_bar)  
          }
          if(f==12){
            for(i in seq(11,61,10)){
              points(x=rep(frequency[i],2),y=c(spectrum[i],coord),type="l",col=col_seasonal,lwd=lwd_seasonal)  
            }
            
            for(i in c(43,53)){
              points(x=rep(frequency[i],2),y=c(spectrum[i],coord),type="l",col=col_td,lwd=lwd_td)  
            }
            
          }else if(f==4){
            for(i in c(31,61)){
              points(x=rep(frequency[i],2),y=c(spectrum[i],coord),type="l",col=col_seasonal,lwd=lwd_seasonal)  
            }
          }
        }else{
          abline(v=(1:(f/2))*1/f,col=col_seasonal,lwd=lwd_seasonal)
          if(f==12)
            abline(v=frequency[c(43,53)],col=col_td,lwd=lwd_td)
          coord <- par("usr")[3]
          for(i in 1:length(frequency)){
            points(x=rep(frequency[i],2),y=c(spectrum[i],coord),type="l",col=col_bar,lwd=lwd_bar)  
          }
        }
        if(f==12){
          aT <- frequency[c(seq(11,61,10),43,53)]
          aL <- paste(round(aT*12),"/",12,sep="")
          aL[(length(aL)-1):length(aL)] <- c("0.348","0.432")
          axis(side=1,at=aT,labels=aL,las=2)
        }else{
          aT <- frequency[c(31,61)]
          aL <- paste(aT*12,"/",12,sep="")
          axis(side=1,at=aT,labels=aL)
        }
        par(mar = c(0, 0, 0, 0)) 
        plot(frequency,spectrum,type = "n", axes = FALSE, ann = FALSE,...)
        if(plot_legend){
          if(f==12)
            legend("center",legend=c("Spectrum","Seasonal Freq.","Trading Day Freq."),
                lty=rep(1,3),col=c(col_bar,col_seasonal,col_td),bg="white",horiz=legend_horiz,bty=legend_bty,lwd=lwd_bar)
          else
            legend("center",legend=c("Spectrum","Seasonal Freq."),lty=rep(1,2),
                col=c(col_bar,col_seasonal),bg="white",horiz=legend_horiz,bty=legend_bty,lwd=lwd_bar)
        }
        gp <- gp[-which(names(gp)=="page")]
        par(gp)
      },finally=par(gp))
}

plot_rsd_acf <- function(out,which="acf",xlab="Lag",ylab="ACF",
    main="default",
    col_acf="black",col_ci="blue",lt_ci=2,ylim="default"){
  x <-out$dg
  if(which=="acf")
    which <- "rsd.acf"
  else if(which=="pacf")
    which <- "rsd.pacf"
  else if(which=="acf2")
    which <- "rsd.acf2"		
  #lwd_bar=4,plot_legend=TRUE){
  if(main=="default"){
    if(which=="rsd.acf"){main <- "Autocorrelations of the Residuals"}
    else if(which=="rsd.pacf"){main <- "Partial Autocorrelations of the Residuals"}        
    else if(which=="rsd.acf2"){main <- "Autocorrelations of the Squared Residuals"}
  }
  if(ylim=="default"){
    ylim<-c(-max(x[[which]][[grep("sample",names(x[[which]]),value=TRUE)]],2*x[[which]][[grep("stderr",names(x[[which]]),value=TRUE)]]),max(x[[which]][[grep("sample",names(x[[which]]),value=TRUE)]],2*x[[which]][[grep("stderr",names(x[[which]]),value=TRUE)]]))
  }
  if(which=="rsd.pacf")
    ylab="Partial ACF"
  plot(x[[which]]$lag,x[[which]][[grep("sample",names(x[[which]]),value=TRUE)]],type="h",xlab=xlab,ylab=ylab,main=main,col=col_acf,ylim=ylim,xaxt="n")
  abline(h=0,col=col_acf)
  if(length(x[[which]]$lag)%in%c(12,24)){
    aT <- c(6,12,18,24)
    axis(side=1,at=aT)  
  }else{
    aT <- c(4,8,12,16)
    axis(side=1,at=aT)
  }
  lines(x[[which]]$lag,2*x[[which]][[grep("stderr",names(x[[which]]),value=TRUE)]],type="l",col=col_ci,lty=lt_ci)
  lines(x[[which]]$lag,-2*x[[which]][[grep("stderr",names(x[[which]]),value=TRUE)]],type="l",col=col_ci,lty=lt_ci)
}

