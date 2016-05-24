genplot<-function(x,y=NA,input="openair",pollutant=NA,plot=TRUE,distr="lnorm",xlim=c(NA,NA),ylim=c(NA,NA),n=NA,col="#A52375",col.points="black",pch=1,il="ls",ci=TRUE,r=0.95,columns=2,xlab="Date",ylab="Concentration",main=NA) {
  
  ## Selection of relevant compound(s).
  # If x is defined as genasis-type data frame, load the compounds.
  if (class(x)=="data.frame"&input=="genasis") {
    compounds<-unique(x[,2])
  }
  
  # If x is defined as openair-type data frame, load the compounds.
  if (class(x)=="data.frame"&input=="openair") {
    compounds<-colnames(x)[-which(is.element(colnames(x),c("date","date_end","temp","wind","note")))]
  }
  
  # If x is defined as numeric vector, set compounds=NA.
  if (class(x)!="data.frame") {
    compounds<-NA
    if(!is.na(pollutant)&pollutant!="") {
      compounds<-pollutant[1]
    }
  }
  
  if (length(pollutant)==1) {
    if(is.na(pollutant)|pollutant=="") {
    pollutant<-compounds
    }
  }
  comp<-pollutant[which(is.element(pollutant,compounds))] # Compounds, which will be used.
  if (length(pollutant)>length(comp)) {
    warning(paste0("One or more pollutants (",pollutant[which(!is.element(pollutant,comp))],") was not recognized."))
  }
  
  # Settings of plotting output.
  if (plot==TRUE&length(comp)>1) {
    par(mfrow=c(ceiling(length(comp)/columns),columns))
  }
  if (plot==TRUE&length(comp)<=1) {
    par(mfrow=c(1,1))
  }
  if (length(comp)>6) {
    warning(paste0("There is a lot of compounds to plot in the dataset x (",length(comp),"). The plot could be crowded on some devices."))
  }
  
  # Output variables.
  if (!is.na(n)) {
    m<-n
  } else {
    if (class(x)=="data.frame") {m<-nrow(x)}
    if (class(x)!="data.frame") {m<-length(x)}
  }
  gslope<-c()
  gintercept<-c()
  gbelt.x<-as.data.frame(matrix(NA,0,m))
  gline<-as.data.frame(matrix(NA,0,m))
  glower<-as.data.frame(matrix(NA,0,m))
  gupper<-as.data.frame(matrix(NA,0,m))
  
   
  ## Cycle over compound(s).
  for (compound in comp) {
    
    ## Preparation of the data for the plot.
    # If x is defined as genasis-type data frame, assignes the variables.
    if (class(x)=="data.frame"&input=="genasis") { 
      valu<-as.numeric(x[which(x[,2]==compound),1])
      date_start<-x[which(x[,2]==compound),3]
      date_end<-x[which(x[,2]==compound),4]
    }
    
    # If x is defined as openair-type data frame, assignes the variables.
    if (class(x)=="data.frame"&input=="openair") {
      if (!is.element("date",colnames(x))){
        stop("Column \"date\" not found.")
      } 
      valu<-as.numeric(x[,compound])
      date_start<-x[,"date"]
      if (is.element("date_end",colnames(x))) {
        date_end<-x[,"date_end"]
      } else {
        date_end<-x[,"date"]
      }
    }
    
    # If x is defined as numeric vector, assignes the variables.
    if (class(x)!="data.frame") {
      if (length(x)!=length(y)) {
        stop(paste0("The length of vectors of concentrations and dates differ."))
      }
      valu<-as.numeric(x)
      date_start<-y
      date_end<-y
    }
    
    
    # If date is given only as the year, shifts the value to 1st june.
    if (max(nchar(as.character(date_start)),na.rm=TRUE)==4) {
      date_start<-paste0(x,"-06-01")
    }
    if (max(nchar(as.character(date_end)),na.rm=TRUE)==4) {
      date_end<-paste0(x,"-06-01")
    }
    
    # Main headers.    
    if (is.na(main)) {
      plotmain<-as.character(compound)
    } else {
      if (class(x)=="data.frame") {
        plotmain<-paste0(main," (",as.character(compound),")")
      } else
        plotmain<-main     
    }
    
    # Transforms concentration values, if option distr=="lnorm"
    if (distr=="lnorm") {
      valu<-log(valu)
    }
    
    # Combines date_start and date_end to one date.
    date_start<-as.Date(as.character(date_start))
    date_end  <-as.Date(as.character(date_end))
    date      <-as.Date((as.numeric(date_start)+as.numeric(date_end))/2,origin="1970-01-01")
    
    valid<-which(!is.na(valu)&!is.na(date))
    valu<-valu[valid]
    date<-date[valid]
    
    
    if (!is.na(max(xlim))) {
      xlim<-as.numeric(as.POSIXct(xlim, format="%Y-%m-%d",tz="GMT"))
    }
    ylim<-as.numeric(ylim)
    
    x0<-as.numeric(as.POSIXct(date, format="%Y-%m-%d",tz="GMT"))    
    y<-valu
    
    xx<-(x0-mean(x0))/(365.256363051*86400)
    xlim<-(xlim-mean(x0))/(365.256363051*86400)
    
    # Least squares interlacing line.
    if (il=="ls") {
      model<-lm(y~xx)
      slope<-as.numeric(unlist(model[1])[2])
      intercept<-as.numeric(unlist(model[1])[1])
    }    
    
    # Theil-Sen interlacing line.
    if (il=="ts") {
      sl<-rep(NA,(length(xx)^2-length(xx))/2)
      ic<-rep(NA,(length(xx)^2-length(xx))/2)
      for (i in 1:(length(xx)-1)) {
        for (j in (i+1):length(xx)) {
          sl[j+length(xx)*(i-1)-((i^2+i)/2)]<-(y[j]-y[i])/(xx[j]-xx[i])
          ic[j+length(xx)*(i-1)-((i^2+i)/2)]<-(xx[j]*y[i]-xx[i]*y[j])/(xx[j]-xx[i])
        }
      }
      slope<-median(sl,na.rm=TRUE)
      intercept<-median(ic,na.rm=TRUE)
    }
    
    
    # Definition of time points
    if (is.na(n)) {
      belt=xx
      if (!is.na(xlim[1])&xlim[1]<min(xx)) {
        breaks<-ceiling((min(xx)-xlim[1])/((max(xx)-min(xx))/(length(xx)-1)))
        addleft<-seq(from=xlim[1],to=min(xx),length.out=breaks+1)
        belt<-c(addleft[-length(addleft)],belt)      
      }
      if (!is.na(xlim[2])&xlim[2]>max(xx)) {
        breaks<-ceiling((xlim[2]-max(xx))/((max(xx)-min(xx))/(length(xx)-1)))
        addright<-seq(from=max(xx),to=xlim[2],length.out=breaks+1)
        belt<-c(belt,addright[-1])      
      }
    }
    
    
    first<-min(xx)-0.001*(max(xx)-min(xx)) # Adds 0.1% reserve to pretend problems with rounding inside min & max functions.
    last<-max(xx)+0.001*(max(xx)-min(xx)) # Adds 0.1% reserve to pretend problems with rounding inside min & max functions.
    
    if (!is.na(xlim[1])) {
      first<-xlim[1]
    }
    if (!is.na(xlim[2])) {
      last<-xlim[2]
    }
    if (!is.na(n)) {
      belt<-(last-first)/(n-1)*rep(0:(n-1))+first
    }
    
    belt<-belt[which(belt>=first&belt<=last)]
    belt.x<-as.Date((belt+mean(xx))*(365.256363051)+mean(x0)/86400,origin="1970-01-01")
    variance<-sqrt(sum((y-(intercept+slope*xx))^2)/(length(xx)-2)*(1/length(xx)+(belt-mean(xx))^2/sum((xx-mean(xx))^2)))
    
    # Least squares interlacing line.
    if (il=="ls") {
      line<-intercept+slope*belt
      upper<-intercept+slope*belt+qt(1-(1-r)/2,length(xx))*variance
      lower<-intercept+slope*belt-qt(1-(1-r)/2,length(xx))*variance
      grmin<-intercept+slope*belt-qt(1-(1-0.95)/2,length(xx))*variance
      grmax<-intercept+slope*belt+qt(1-(1-0.95)/2,length(xx))*variance
    }
    
    # Theil-Sen interlacing line.
    if (il=="ts") {
      line<-intercept+slope*belt
      upper<-pmax(quantile(ic,(1+r)/2)+quantile(sl,(1-r)/2)*belt,quantile(ic,(1+r)/2)+quantile(sl,(1+r)/2)*belt)
      lower<-pmin(quantile(ic,(1-r)/2)+quantile(sl,(1-r)/2)*belt,quantile(ic,(1-r)/2)+quantile(sl,(1+r)/2)*belt)
      grmin<-pmin(quantile(ic,0.975)+quantile(sl,0.975)*belt,quantile(ic,0.025)+quantile(sl,0.025)*belt)
      grmax<-pmax(quantile(ic,0.975)+quantile(sl,0.975)*belt,quantile(ic,0.025)+quantile(sl,0.025)*belt)
    }
    

    # Reverse transformation.    
    if (distr=="lnorm") {
      line<-exp(line)
      upper<-exp(upper)
      lower<-exp(lower)
      grmax<-exp(grmax)
      grmin<-exp(grmin)
      y<-exp(y)
    } 
    
    # Plot.
    if (length(valu)==1){
      if (is.na(ylim[1])) {
        ylim[1]=c(0.9*valu[1])
      }
      if (is.na(ylim[2])){
        ylim[2]=c(1.1*valu[1])
      }
    }
    
    plotylim<-ylim
    if (plot==TRUE) {
      if (is.na(max(ylim))) {
        if(ci==TRUE&length(valu)>2) {
          plotylim=c(min(lower,y),max(upper,y))
        }
        if(ci=="gradient"&length(valu)>2) {
          plotylim=c(min(grmin,y),max(grmax,y))
        }
        if ((ci!=TRUE&ci!="gradient")|length(valu)<=2) {
          plotylim=c(min(line,y),max(line,y))
        }
      }
      
      plot(belt.x,line,cex=0,xlab=xlab,ylab=ylab,ylim=plotylim,main=plotmain,pch=pch)
      
      # Plot interlacing curve.
      if (il!="none"&il!=""&!is.na(il)) {
        lines(belt.x,line,col=col,lwd=3)
      }
      
      # Plot confidence interval.
      if (ci==TRUE&length(valu)>2) {
        for (i in 2:length(belt)) {
          polygon(c(belt.x[i-1],belt.x[i-1],belt.x[i],belt.x[i]),c(upper[i-1],lower[i-1],lower[i],upper[i]),border=NA,col=paste0(rgb(t(col2rgb(col)/255)),"80"))
        }
      }
      
      if (ci=="gradient"&length(valu)>2) {
        for (g in 1:99) {
          rg<-1-0.01*g
           
          if (il=="ls") {
            if (distr=="lnorm") {
              line<-exp(intercept+slope*belt)
              upper<-exp(intercept+slope*belt+qt(1-(1-rg)/2,length(xx))*variance)
              lower<-exp(intercept+slope*belt-qt(1-(1-rg)/2,length(xx))*variance)
            } else {
              line<-intercept+slope*belt
              upper<-intercept+slope*belt+qt(1-(1-rg)/2,length(xx))*variance
              lower<-intercept+slope*belt-qt(1-(1-rg)/2,length(xx))*variance
            }
          }
          if (il=="ts") {
            if (distr=="lnorm") {
              line<-exp(intercept+slope*belt)
              upper<-exp(pmax(quantile(ic,(1+rg)/2)+quantile(sl,(1-rg)/2)*belt,quantile(ic,(1+rg)/2)+quantile(sl,(1+rg)/2)*belt))
              lower<-exp(pmin(quantile(ic,(1-rg)/2)+quantile(sl,(1-rg)/2)*belt,quantile(ic,(1-rg)/2)+quantile(sl,(1+rg)/2)*belt))
            } else {
              line<-intercept+slope*belt
              upper<-pmax(quantile(ic,(1+rg)/2)+quantile(sl,(1-rg)/2)*belt,quantile(ic,(1+rg)/2)+quantile(sl,(1+rg)/2)*belt)
              lower<-pmin(quantile(ic,(1-rg)/2)+quantile(sl,(1-rg)/2)*belt,quantile(ic,(1-rg)/2)+quantile(sl,(1+rg)/2)*belt)
            }
          }
          
          for (i in 2:length(belt)) {
            polygon(c(belt.x[i-1],belt.x[i-1],belt.x[i],belt.x[i]),c(upper[i-1],lower[i-1],lower[i],upper[i]),border=NA,col=paste0(rgb(t(col2rgb(col)/255)),"06"))
          }
        }
      }
      points(as.Date(date,origin="1970-01-01"),y,col=col.points,pch=pch) 
    }
    gslope<-c(gslope,slope/(365.256363051*86400))
    gintercept<-c(gintercept,intercept)
    gbelt.x<-rbind(gbelt.x,belt.x)
    gline<-rbind(gline,line)
    glower<-rbind(glower,lower)
    gupper<-rbind(gupper,upper)
  }
par(mfrow=c(1,1))
invisible(list(slope=gslope,intercept=gintercept,belt=gbelt.x,line=gline,lower=glower,upper=gupper))
}