genanaggr<-function(x,y=NA,input="openair",output=NA,pollutant=NA,method="mean",minn=4,gap=3,show.flagged=FALSE) {
  gmean<-function(x) {
    geom<-exp(mean(log(x),na.rm=TRUE))
    return(geom)
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
  
  
  # Output type selection.
  if (is.na(output)&class(x)=="data.frame") {
    output<-input
  }
  
  
  ## Creating column variables.
  # If x is defined as genasis-type data frame, set column vectors.
  if (class(x)=="data.frame"&input=="genasis") {
    valu      <-as.numeric(x[,1])
    comp      <-as.character(x[,2])
    date_start<-gendate(x[,3])
    date_end  <-gendate(x[,4])
    temp      <-as.numeric(x[,5])
    wind      <-as.numeric(x[,6])
    
  }
  
  # If x is defined as openair-type data frame, set column vectors.  
  if (class(x)=="data.frame"&input=="openair") {
    valu      <-c()
    comp      <-c()
    date_start<-c()
    date_end  <-c()
    temp      <-c()
    wind      <-c()
    for (compound in compos) {
      valu      <-c(valu,as.numeric(x[,compound]))
      comp      <-c(comp,as.character(rep(compound,nrow(x))))
      date_start<-as.Date(c(as.character(date_start),as.character(x[,"date"])))
      if (is.element("date_end",colnames(x))) {
        date_end<-as.Date(c(as.character(date_end),as.character(x[,"date_end"])))} else {
          date_end<-as.Date(c(as.character(date_end),as.character(x[,"date"])))
        }
      if (is.element("temp",colnames(x))) {
        temp<-c(temp,as.numeric(x[,"temp"]))} else {
          temp<-c(temp,rep(NA,nrow(x)))
        }
      if (is.element("wind",colnames(x))) {
        wind<-c(wind,as.numeric(x[,"wind"]))} else {
          wind<-c(wind,rep(NA,nrow(x)))
        }
    }
  }
  
  # If x is defined as numeric vector, set column vectors.
  if (class(x)!="data.frame") {
    valu      <-x
    comp      <-rep(pollutant,length(x))
    date_start<-gendate(y)
    date_end  <-gendate(y)
    temp      <-rep(NA,length(x))
    wind      <-rep(NA,length(x)) 
  }
  
  
  ## Combines date_start and date_end to one date.
  date<-as.Date((as.numeric(date_start)+as.numeric(date_end))/2,origin="1970-01-01")
  
 
  ## Aggregation itself.
  rvalu<-c()
  rcomp<-c()
  rdate_start<-c()
  rdate_end<-c()
  rtemp<-c()
  rwind<-c()
  rnote<-c()
  
  for (compound in compos) {
    svalu<-valu[which(comp==compound)]
    scomp<-comp[which(comp==compound)]
    sdate<-date[which(comp==compound)]
    stemp<-temp[which(comp==compound)]
    swind<-wind[which(comp==compound)]
    
    if (class(x)!="data.frame") {
      svalu<-x
      scomp<-pollutant
      sdate<-y
      stemp<-NA
      swind<-NA
    }
    
    valid<-which(!is.na(svalu)&!is.na(sdate))
    svalu<-svalu[valid]
    scomp<-scomp[valid]
    sdate<-sdate[valid]
    stemp<-stemp[valid]
    swind<-swind[valid]
    
    for (j in initial:final) {
      ssvalu<-svalu[which(as.numeric(substr(sdate,1,4))==j)]
      sscomp<-scomp[which(as.numeric(substr(sdate,1,4))==j)]
      ssdate<-sdate[which(as.numeric(substr(sdate,1,4))==j)]
      sstemp<-stemp[which(as.numeric(substr(sdate,1,4))==j)]
      sswind<-swind[which(as.numeric(substr(sdate,1,4))==j)]
      
      ssdate<-ssdate[order(ssdate)]
      sstemp<-sstemp[order(ssdate)]
      sswind<-sswind[order(ssdate)]
      
      intervals0<-c(as.Date(paste0(j,"-01-01")),as.Date(ssdate),as.Date(paste0(j,"-12-31")))
      intervals<-as.numeric(intervals0[-1]-intervals0[-length(intervals0)])
      
      report<-NA
      if (length(ssdate[!is.na(ssdate)])<minn) {
        report<-paste0("Less than ",minn," values in the year.")
      }
      if (max(intervals,na.rm=TRUE)>gap*365/length(intervals)) {
        report<-"Measurements are too irregularly distributed in time."
      }
      if (length(ssdate)<1) {
        report<-"No valid value within the year."
      }
      
      if (!is.na(report)&show.flagged==FALSE) {
        rvalu<-c(rvalu,NA)
      } else {
        rvalu<-c(rvalu,unlist(lapply(list(ssvalu),FUN=method)))
      }
      rcomp<-c(rcomp,compound)
      rdate_start<-c(rdate_start,intervals0[1])
      rdate_end<-c(rdate_end,intervals0[length(intervals0)])
      rtemp<-c(rtemp,mean(sstemp,na.rm=TRUE))
      rwind<-c(rwind,mean(sswind,na.rm=TRUE))
      rnote<-c(rnote,report)
    }
  }
  
  ## Generation of desired data type.
  # Output in genasis mode.
  if (!is.na(output)&output=="genasis") {
    res<-data.frame(rvalu,rcomp,as.Date(rdate_start,origin="1970-01-01"),as.Date(rdate_end,origin="1970-01-01"),rtemp,rwind,rnote)
    colnames(res)<-c("valu","comp","date_start","date_end","temp","wind","note")
  }
  
  # Output in openair mode.
  if (!is.na(output)&output=="openair") {
    unidates<-as.data.frame(matrix(NA,0,2))
    for (a in unique(rdate_start)) {
      for (b in unique(rdate_end[which(rdate_start==a)])) {
        unidates<-rbind(unidates,as.Date(c(a,b),origin="1970-01-01"))
      }
    }
    
    unidates[,1]<-as.Date(unidates[,1],origin="1970-01-01")
    unidates[,2]<-as.Date(unidates[,2],origin="1970-01-01")
    
    res<-as.data.frame(matrix(NA,nrow(unidates),length(compos)+5))
    colnames(res)<-c("date","date_end","temp","wind","note",compos)
    
    
    for (i in 1:nrow(unidates)) {
      res[i,1]<-unidates[i,1]
      res[i,2]<-unidates[i,2]
      res[i,3]<-rtemp[which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2]&rcomp==compound)][1]
      res[i,4]<-rwind[which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2]&rcomp==compound)][1]     
      if (paste(rcomp[which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2]&!is.na(rnote))], collapse=', ')=="") {
        res[i,5]<-NA
      } else {
        res[i,5]<-paste0("Problems with: ",paste(rcomp[which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2]&!is.na(rnote))], collapse=', '),".")
      }
      for (compound in compos) {
        res[i,compound]<-rvalu[which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2]&rcomp==compound)][1]
      }
      if (length(which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2]&rcomp==compound))>1) {
        warning(paste0("There was more than 1 record with the same start and end date for ",compound,", only the first record was used."))
        
      }
    }
    res[,1]<-as.Date(res[,1],origin="1970-01-01")
    res[,2]<-as.Date(res[,2],origin="1970-01-01")
  }
  
  # Output in vector-mode.
  if (class(x)!="data.frame"&is.na(output)) {
    res<-rvalu 
  }
  
  return(res)
}