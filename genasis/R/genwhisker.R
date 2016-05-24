genwhisker<-function(x,y=NA,input="openair",method="mqm",pollutant=NA,distr="norm",by.years=FALSE,col="#A52375",legend=TRUE,xlab="",ylab="Concentration",main=NA) {
  quantile05<-function(x) {
    quantile(x, 0.05, names=FALSE, na.rm=TRUE)
  }
  quantile25<-function(x) {
    quantile(x, 0.25, names=FALSE, na.rm=TRUE)
  }
  quantile75<-function(x) {
    quantile(x, 0.75, names=FALSE, na.rm=TRUE)
  }
  quantile95<-function(x) {
    quantile(x, 0.95, names=FALSE, na.rm=TRUE)
  }
  trimean<-function(x) {
    (quantile25+2*median(x,na.rm=TRUE)+quantile75)/4
  }
  
  ## Selection of relevant compound(s).
  # If x is defined as genasis-type data frame, load the compounds.
  if (class(x)=="data.frame"&input=="genasis") {
    compounds<-as.character(unique(x[,2]))
    initial<-as.numeric(substr(min(gendate(x[,3]),na.rm=TRUE)+(min(gendate(x[,4]),na.rm=TRUE)-min(gendate(x[,3]),na.rm=TRUE))/2,1,4))
    final<-as.numeric(substr(max(gendate(x[,4]),na.rm=TRUE)-(max(gendate(x[,4]),na.rm=TRUE)-max(gendate(x[,3]),na.rm=TRUE))/2,1,4))
  }
  
  # If x is defined as openair-type data frame, load the compounds.
  if (class(x)=="data.frame"&input=="openair") {
    compounds<-colnames(x)[-which(is.element(colnames(x),c("date","date_end","temp","wind","note")))]
    x[,"date"]<-gendate(x[,"date"])
    if (is.element("date_end",colnames(x))) {
      x[,"date_end"]<-gendate(x[,"date_end"])
      initial<-as.numeric(substr(min(x[,"date"],na.rm=TRUE)+(min(x[,"date_end"],na.rm=TRUE)-min(x[,"date"],na.rm=TRUE))/2,1,4))
      final<-as.numeric(substr(max(x[,"date_end"],na.rm=TRUE)-(max(x[,"date_end"],na.rm=TRUE)-max(x[,"date"],na.rm=TRUE))/2,1,4))
    } else {
      initial<-min(as.numeric(substr(x[,"date"],1,4)),na.rm=TRUE)
      final<-max(as.numeric(substr(x[,"date"],1,4)),na.rm=TRUE)
    }
  }
  
  # If x is vector, load the compounds.
  if (class(x)!="data.frame") {
    compounds<-NA
    initial<-min(as.numeric(substr(gendate(y),1,4)),na.rm=TRUE)
    final<-max(as.numeric(substr(gendate(y),1,4)),na.rm=TRUE)
  }
  
  ## Selection of compounds into variable compos.
  if (length(pollutant)==1) {
    if(is.na(pollutant)|pollutant=="") {
      pollutant<-compounds
    }
  }
  if (class(x)=="data.frame") {
    compos<-pollutant[which(is.element(pollutant,compounds))] # Compounds, which will be used.
    
    if (length(pollutant)>length(compos)) {
      warning(paste0("One or more pollutants (",pollutant[which(!is.element(pollutant,compos))],") was not recognized."))
    }
  } else {
    compos<-pollutant
  }
  
  ## Creating column variables.
  # If x is defined as genasis-type data frame, set column vectors.
  if (class(x)=="data.frame"&input=="genasis") {
    valu      <-as.numeric(x[,1])
    comp      <-as.character(x[,2])
    date_start<-gendate(x[,3])
    date_end  <-gendate(x[,4])
  }
  
  # If x is defined as openair-type data frame, set column vectors.  
  if (class(x)=="data.frame"&input=="openair") {
    valu      <-c()
    comp      <-c()
    date_start<-c()
    date_end  <-c()
    for (compound in compos) {
      valu      <-c(valu,as.numeric(x[,compound]))
      comp      <-c(comp,as.character(rep(compound,nrow(x))))
      date_start<-as.Date(c(as.character(date_start),as.character(x[,"date"])))
      if (is.element("date_end",colnames(x))) {
        date_end<-as.Date(c(as.character(date_end),as.character(x[,"date_end"])))} else {
          date_end<-as.Date(c(as.character(date_end),as.character(x[,"date"])))
        }
    }
  }
  
  # If x is defined as numeric vector, set column vectors.
  if (class(x)!="data.frame") {
    valu      <-x
    comp      <-rep(pollutant,length(x))
    date_start<-gendate(y)
    date_end  <-gendate(y)
  }
  
  
  ## Combines date_start and date_end to one date.
  date<-as.Date((as.numeric(date_start)+as.numeric(date_end))/2,origin="1970-01-01")
  
  # Transforms concentration values, if option distr=="lnorm"
  if (distr=="lnorm") {
    valu<-log(valu)
  }
  
  ## Plotting.
  highest<-max(valu,na.rm=TRUE)
  lowest<-min(valu,na.rm=TRUE)
  unit<-(highest-lowest)/100
  
  if(method=="mqm") {
    f1<-"median"
    f2<-"quantile25"
    f3<-"quantile75"
    f4<-"min"
    f5<-"max"
  }
  if(method=="tqm") {
    f1<-"trimean"
    f2<-"quantile25"
    f3<-"quantile75"
    f4<-"min"
    f5<-"max"
  }
  if(method=="mqq") {
    f1<-"median"
    f2<-"quantile25"
    f3<-"quantile75"
    f4<-"quantile05"
    f5<-"quantile95"
  }
  
  ## Plot
  par(mar=c(9,4,4,2),mfrow=c(1,1))
  if (distr=="lnorm") {
    if (by.years==TRUE&legend==TRUE) {
      plot(c(1:(length(compos)+1)),c(rep(lowest,length(compos)),highest+0.2*(highest-lowest)),cex=0,xaxt="n",yaxt="n",xlab=xlab,ylab=ylab,main=main)     
    } else {
      plot(c(1:(length(compos)+1)),c(rep(lowest,length(compos)),highest),cex=0,xaxt="n",yaxt="n",xlab=xlab,ylab=ylab,main=main)
    }
    axis(2,at=axis(2,labels=NA),round(exp(axis(2,labels=NA)),2))
  } else {
    if (by.years==TRUE&legend==TRUE) {
      plot(c(1:(length(compos)+1)),c(rep(lowest,length(compos)),highest+0.2*(highest-lowest)),cex=0,xaxt="n",xlab=xlab,ylab=ylab,main=main)
    } else {
    plot(c(1:(length(compos)+1)),c(rep(lowest,length(compos)),highest),cex=0,xaxt="n",xlab=xlab,ylab=ylab,main=main)
    }
  }
  axis(1,at=c(1:length(compos)+0.5),labels=compos,las=2)
  
  for (compound in compos) {
    
    if (class(x)=="data.frame") {
      i<-which(compos==compound)
    } else {
      i<-1
    }
    
    svalu<-valu[which(comp==compound)]
    scomp<-comp[which(comp==compound)]
    sdate<-date[which(comp==compound)]
    
    if (class(x)!="data.frame") {
      svalu<-valu
      scomp<-pollutant
      sdate<-y
    }
    
    valid<-which(!is.na(svalu)&!is.na(sdate))
    svalu<-svalu[valid]
    scomp<-scomp[valid]
    sdate<-sdate[valid]
    
    # One-coloured plot.
    if (by.years==FALSE) {
      # Central point.
      rect(i-0.4+0.5,apply(as.matrix(svalu),2,FUN=f1)-0.5*unit,i+0.4+0.5,apply(as.matrix(svalu),2,FUN=f1)+0.5*unit,border=NA,col=paste0(rgb(t(col2rgb(col)/255)),"84"))
      # Upper vertical bar.
      rect(i-0.1+0.5,apply(as.matrix(svalu),2,FUN=f3),i+0.1+0.5,apply(as.matrix(svalu),2,FUN=f5)-0.5*unit,border=NA,col=paste0(rgb(t(col2rgb(col)/255)),"60"))
      # Lower vertical bar.
      rect(i-0.1+0.5,apply(as.matrix(svalu),2,FUN=f4)+0.5*unit,i+0.1+0.5,apply(as.matrix(svalu),2,FUN=f2),border=NA,col=paste0(rgb(t(col2rgb(col)/255)),"60"))
      # Box.
      rect(i-0.4+0.5,apply(as.matrix(svalu),2,FUN=f2),i+0.4+0.5,apply(as.matrix(svalu),2,FUN=f3),border=NA,col=paste0(rgb(t(col2rgb(col)/255)),"60"))
      # Lower whisker.
      rect(i-0.4+0.5,apply(as.matrix(svalu),2,FUN=f4)-0.5*unit,i+0.4+0.5,apply(as.matrix(svalu),2,FUN=f4)+0.5*unit,border=NA,col=paste0(rgb(t(col2rgb(col)/255)),"60"))
      # Upper whisker.
      rect(i-0.4+0.5,apply(as.matrix(svalu),2,FUN=f5)-0.5*unit,i+0.4+0.5,apply(as.matrix(svalu),2,FUN=f5)+0.5*unit,border=NA,col=paste0(rgb(t(col2rgb(col)/255)),"60"))
    }
    
    # Multi-coloured plot.
    
    if (by.years==TRUE) {
      
      # Colours generation.
      if (length(col)!=final-initial+1) {
        col<-hsv((1:(final-initial+1))/(final-initial+1),1,1)
      }
      
      #
      for (j in initial:final) {
        
        ssvalu<-svalu[which(as.numeric(substr(sdate,1,4))==j)]
        
        if (length(ssvalu)>0) {
          # Central point.
          rect(i-0.4+0.5,apply(as.matrix(ssvalu),2,FUN=f1)-0.5*unit,i+0.4+0.5,apply(as.matrix(ssvalu),2,FUN=f1)+0.5*unit,border=NA,col=paste0(rgb(t(col2rgb(col[j-initial+1])/255)),"84"))
          # Upper vertical bar.
          rect(i-0.1+0.5,apply(as.matrix(ssvalu),2,FUN=f3),i+0.1+0.5,apply(as.matrix(ssvalu),2,FUN=f5)-0.5*unit,border=NA,col=paste0(rgb(t(col2rgb(col[j-initial+1])/255)),"60"))
          # Lower vertical bar.
          rect(i-0.1+0.5,apply(as.matrix(ssvalu),2,FUN=f4)+0.5*unit,i+0.1+0.5,apply(as.matrix(ssvalu),2,FUN=f2),border=NA,col=paste0(rgb(t(col2rgb(col[j-initial+1])/255)),"60"))
          # Box.
          rect(i-0.4+0.5,apply(as.matrix(ssvalu),2,FUN=f2),i+0.4+0.5,apply(as.matrix(ssvalu),2,FUN=f3),border=NA,col=paste0(rgb(t(col2rgb(col[j-initial+1])/255)),"60"))
          # Lower whisker.
          rect(i-0.4+0.5,apply(as.matrix(ssvalu),2,FUN=f4)-0.5*unit,i+0.4+0.5,apply(as.matrix(ssvalu),2,FUN=f4)+0.5*unit,border=NA,col=paste0(rgb(t(col2rgb(col[j-initial+1])/255)),"60"))
          # Upper whisker.
          rect(i-0.4+0.5,apply(as.matrix(ssvalu),2,FUN=f5)-0.5*unit,i+0.4+0.5,apply(as.matrix(ssvalu),2,FUN=f5)+0.5*unit,border=NA,col=paste0(rgb(t(col2rgb(col[j-initial+1])/255)),"60"))
        }
      }
      # Legend
      if (legend==TRUE) {
        for (j in (initial:final-initial+1)) {
          legend((j-1)/(final-initial+1)*length(compos)*0.9+1,highest+0.3*(highest-lowest),(initial:final)[j],horiz=TRUE,text.col=paste0(rgb(t(col2rgb(col[j])/255)),"60"),bty="n")
        }
      }
    }
  }
  par(mar=c(5,4,4,2))
}