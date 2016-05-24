genstatistic<-function(x,y=NA,input="openair",pollutant=NA,by.years=FALSE) {
  gmean<-function(x) {
    geom<-exp(mean(log(x[which(!is.na(x)&x>0)]),na.rm=TRUE))
    if (min(x,na.rm=TRUE)<=0){
      geom<-NA
    }
    return(geom)
  }
  gsd<-function(x) {
    geom<-exp(sd(log(x[which(!is.na(x)&x>0)]),na.rm=TRUE))
    if (min(x,na.rm=TRUE)<=0){
      geom<-NA
    }
    return(geom)
  }
  senthail<-function(x,y) {
    x<-as.numeric(as.POSIXct(x, format="%Y-%m-%d",tz="GMT"))/86400
    sl<-rep(NA,(length(x)^2-length(x))/2)
    ic<-rep(NA,(length(x)^2-length(x))/2)
    for (i in 1:(length(x)-1)) {
      for (j in (i+1):length(x)) {
        sl[j+length(x)*(i-1)-((i^2+i)/2)]<-(y[j]-y[i])/(x[j]-x[i])
        ic[j+length(x)*(i-1)-((i^2+i)/2)]<-(x[j]*y[i]-x[i]*y[j])/(x[j]-x[i])
      }
    }
    slope<-median(sl,na.rm=TRUE)
    intercept<-median(ic,na.rm=TRUE)
    return(list(slope=slope,intercept=intercept))
  }
  pearcor<-function(x,y){
    if (length(x[which(!is.na(x)&!is.na(y))])>2) {
      cc<-suppressWarnings(cor.test(x,as.numeric(gendate(y)),method="pearson",alternative="two.sided")$estimate)
      p<-suppressWarnings(cor.test(x,as.numeric(gendate(y)),method="pearson",alternative="two.sided")$p.value)
    } else {
      cc<-NA
      p<-NA
    }
    return(list(cc=cc,p=p))
  }
  spearcor<-function(x,y){
    if (length(x[which(!is.na(x)&!is.na(y))])>2) {
      cc<-suppressWarnings(cor.test(x,as.numeric(gendate(y)),method="spearman",alternative="two.sided")$estimate)
      p<-suppressWarnings(cor.test(x,as.numeric(gendate(y)),method="spearman",alternative="two.sided")$p.value)
    } else {
      cc<-NA
      p<-NA
    }
    return(list(cc=cc,p=p))
  }
  mktest<-function(x){
    if (length(x[which(!is.na(x))])>3) {
      cc<-suppressWarnings(MannKendall(x))$tau[1]
      p<-suppressWarnings(MannKendall(x))$sl[1]
    } else {
      cc<-NA
      p<-NA
    }
    return(list(cc=cc,p=p))
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
  
  # If x is defined as numeric vector, set compounds=NA.
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
  
  ## Definition of table of results.
  if (by.years==FALSE) {
    result<-as.data.frame(matrix(NA,length(compos),18))
  }
  if (by.years==TRUE) {
    result<-as.data.frame(matrix(NA,length(compos)*(final-initial+1),19))
  }

  ## Cycle over compound(s).
  for (compound in compos) {
    
    ## Preparation of the data for the plot.
    # If x is defined as genasis-type data frame, assignes the variables.
    if (class(x)=="data.frame"&input=="genasis") { 
      svalu<-as.numeric(x[which(x[,2]==compound),1])
      sdate_start<-gendate(x[which(x[,2]==compound),3])
      sdate_end<-gendate(x[which(x[,2]==compound),4])
    }
    
    # If x is defined as openair-type data frame, assignes the variables.
    if (class(x)=="data.frame"&input=="openair") {
      if (!is.element("date",colnames(x))){
        stop("Column \"date\" not found.")
      } 
      svalu<-as.numeric(x[,compound])
      sdate_start<-gendate(x[,"date"])
      if (is.element("date_end",colnames(x))) {
        sdate_end<-gendate(x[,"date_end"])
      } else {
        sdate_end<-gendate(x[,"date"])
      }
    }
    
    # If x is defined as numeric vector, assignes the variables.
    if (class(x)!="data.frame") {
      if (length(x)!=length(y)) {
        stop(paste0("The length of vectors of concentrations and dates differ."))
      }
      svalu<-x
      sdate_start<-gendate(y)
      sdate_end<-gendate(y)
    }
    
    # Combines date_start and date_end to one date.
    sdate<-as.Date((as.numeric(sdate_start)+as.numeric(sdate_end))/2,origin="1970-01-01")
    
    # Select only records with valid date and value.
    valid<-which(!is.na(svalu)&!is.na(sdate))
    svalu<-svalu[valid]
    sdate<-sdate[valid]
    
    if (by.years==FALSE) {
      i<-which(compos==compound)
      if (is.na(compound)){
        i<-1
      }
      
      result[i,1]<-compound
      result[i,2]<-length(svalu[!is.na(svalu)])
      result[i,3]<-round(mean(svalu,na.rm=TRUE),3)
      result[i,4]<-round(sd(svalu,na.rm=TRUE),3)
      result[i,5]<-round(gmean(svalu),3)
      result[i,6]<-round(gsd(svalu),3)
      result[i,7]<-round(min(svalu,na.rm=TRUE),3)
      result[i,8]<-round(median(svalu,na.rm=TRUE),3)
      result[i,9]<-round(max(svalu,na.rm=TRUE),3)
      result[i,10]<-round(pearcor(svalu,sdate)$cc,3)
      result[i,11]<-round(pearcor(svalu,sdate)$p,3)
      result[i,12]<-round(spearcor(svalu,sdate)$cc,3)
      result[i,13]<-round(spearcor(svalu,sdate)$p,3)
      result[i,14]<-round(mktest(svalu)$cc,3)
      result[i,15]<-round(mktest(svalu)$p,3)
      result[i,16]<-round(unlist(lm(svalu~as.numeric(sdate)+1)[1])[2],3)
      result[i,17]<-round(unlist(senthail(sdate,svalu))[1],3)
      result[i,18]<-round(svalu[length(svalu)]-svalu[1],3)
    }
    if (by.years==TRUE) {
      i<-which(compos==compound)
      if (is.na(compound)){
        i<-1
      }
      for (year in initial:final) {
        j<-year-initial+1
        ssvalu<-svalu[which(as.numeric(substr(sdate,1,4))==year)]
        ssdate<-sdate[which(as.numeric(substr(sdate,1,4))==year)]
        
        result[(i-1)*(final-initial+1)+j,1]<-compound
        result[(i-1)*(final-initial+1)+j,2]<-year
        result[(i-1)*(final-initial+1)+j,3]<-length(ssvalu[!is.na(ssvalu)])
        result[(i-1)*(final-initial+1)+j,4]<-round(mean(ssvalu,na.rm=TRUE),3)
        result[(i-1)*(final-initial+1)+j,5]<-round(sd(ssvalu,na.rm=TRUE),3)
        result[(i-1)*(final-initial+1)+j,6]<-round(gmean(ssvalu),3)
        result[(i-1)*(final-initial+1)+j,7]<-round(gsd(ssvalu),3)
        result[(i-1)*(final-initial+1)+j,8]<-round(min(ssvalu,na.rm=TRUE),3)
        result[(i-1)*(final-initial+1)+j,9]<-round(median(ssvalu,na.rm=TRUE),3)
        result[(i-1)*(final-initial+1)+j,10]<-round(max(ssvalu,na.rm=TRUE),3)
        result[(i-1)*(final-initial+1)+j,10]<-round(pearcor(ssvalu,ssdate)$cc,3)
        result[(i-1)*(final-initial+1)+j,11]<-round(pearcor(ssvalu,ssdate)$p,3)
        result[(i-1)*(final-initial+1)+j,12]<-round(spearcor(ssvalu,ssdate)$cc,3)
        result[(i-1)*(final-initial+1)+j,13]<-round(spearcor(ssvalu,ssdate)$p,3)
        result[(i-1)*(final-initial+1)+j,15]<-round(mktest(ssvalu)$cc,3)
        result[(i-1)*(final-initial+1)+j,16]<-round(mktest(ssvalu)$p,3)
        result[(i-1)*(final-initial+1)+j,17]<-round(unlist(lm(ssvalu~as.numeric(ssdate)+1)[1])[2],3)
        result[(i-1)*(final-initial+1)+j,18]<-round(unlist(senthail(ssdate,ssvalu))[1],3)
        result[(i-1)*(final-initial+1)+j,19]<-round(ssvalu[length(ssvalu)]-ssvalu[1],3)
      }
    }

  }
  
  if (by.years==FALSE) {
    colnames(result)<-c("pollutant","n","mean","sd","geom. mean","geom. sd","min","median","max","Pearson","Pp","Daniels","Dp","Mann-Kendall","MKp","LS slope","TS slope","delta")
  }
  if (by.years==TRUE) {
    colnames(result)<-c("pollutant","year","n","mean","sd","geom. mean","geom. sd","min","median","max","Pearson","Pp","Daniels","Dp","Mann-Kendall","MKp","LS slope","TS slope","delta")
  }
return(list(res=result))  
}