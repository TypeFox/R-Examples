genoutlier<-function(x,y=NA,input="openair",output=NA,method="lm3s",sides=2,pollutant=NA,plot=TRUE,columns=2,col.points="black",pch=1,xlab="Date",ylab="Concentration",main=NA) {
  
  ## Selection of relevant compound(s).
  # If x is defined as genasis-type data frame, load the compounds.
  if (class(x)=="data.frame"&input=="genasis") {
    compounds<-as.character(unique(x[,2]))
  }
  
  # If x is defined as openair-type data frame, load the compounds.
  if (class(x)=="data.frame"&input=="openair") {
    compounds<-colnames(x)[-which(is.element(colnames(x),c("date","date_end","temp","wind","note")))]
  }
  
  # If x is defined as numeric vector, set compounds=pollutant.
  if (class(x)!="data.frame") {
    compounds<-pollutant
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
    date_start<-as.Date(as.character(x[,3]))
    date_end  <-as.Date(as.character(x[,4]))
    if (ncol(x)>4) {
      temp      <-as.numeric(x[,5])
    } else {
      temp<-rep(NA,nrow(x))
    }
    if (ncol(x)>5) {
    wind      <-as.numeric(x[,6])
    } else {
      wind<-rep(NA,nrow(x))
    }
    
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
    date_start<-y
    date_end  <-y
    temp      <-rep(NA,length(x))
    wind      <-rep(NA,length(x)) 
  }
  
  # Settings of plotting output.
  if (plot==TRUE&length(compos)>1) {
    par(mfrow=c(ceiling(length(compounds)/columns),columns))
  }
  if (plot==TRUE&length(compos)<=1) {
    par(mfrow=c(1,1))
  }
  if (plot==TRUE&length(compos)>6) {
    warning(paste0("There is a lot of compounds to plot in the dataset x (",length(compos),"). The plot could be crowded on some devices."))
  }
  
  ## Outliers identification and exclusion
  note<-rep(NA,length(valu))
  glower<-rep(NA,length(compos))
  gupper<-rep(NA,length(compos))

  for (compound in compos) {
    svalu<-valu[which(comp==compound)]
    sdate_start<-date_start[which(comp==compound)]
    sdate_end<-date_end[which(comp==compound)]
    
    if (class(x)!="data.frame") {
      svalu<-x
      sdate_start<-y
      sdate_end<-y
    }
    
    if (method=="m2s") {lower<-mean(svalu,na.rm=TRUE)-2*sd(svalu,na.rm=TRUE);upper<-mean(svalu,na.rm=TRUE)+2*sd(svalu,na.rm=TRUE)}
    if (method=="m3s") {lower<-mean(svalu,na.rm=TRUE)-3*sd(svalu,na.rm=TRUE);upper<-mean(svalu,na.rm=TRUE)+3*sd(svalu,na.rm=TRUE)}
    if (method=="lm2s") {lower<-exp(mean(log(svalu),na.rm=TRUE)-2*sd(log(svalu),na.rm=TRUE));upper<-exp(mean(log(svalu),na.rm=TRUE)+2*sd(log(svalu),na.rm=TRUE))}
    if (method=="lm3s") {lower<-exp(mean(log(svalu),na.rm=TRUE)-3*sd(log(svalu),na.rm=TRUE));upper<-exp(mean(log(svalu),na.rm=TRUE)+3*sd(log(svalu),na.rm=TRUE))}
    if (method=="iqr2") {lower<-(1.5*quantile(svalu,0.25)-0.5*quantile(svalu,0.75));upper<-(1.5*quantile(svalu,0.75)-0.5*quantile(svalu,0.25))}
    if (method=="iqr4") {lower<-(2.5*quantile(svalu,0.25)-1.5*quantile(svalu,0.75));upper<-(2.5*quantile(svalu,0.75)-1.5*quantile(svalu,0.25))}
    if (method=="iqr7") {lower<-(4.0*quantile(svalu,0.25)-3.0*quantile(svalu,0.75));upper<-(4.0*quantile(svalu,0.75)-3.0*quantile(svalu,0.25))}    
       
    for (j in 1:length(svalu)) {
      if (length(svalu)>1) {
        if (sides==1) {
          if ((svalu[j]<=upper)|is.na(svalu[j])) {
            valu[which(comp==compound)][j]<-svalu[j]
            note[which(comp==compound)][j]<-NA
          } else {
            valu[which(comp==compound)][j]<-NA
            if (class(x)!="data.frame") {
              valu[j]<-NA
            }
            if (svalu[j]>upper&!is.na(svalu[j])) {
              note[which(comp==compound)][j]<-paste0("Outlier: upper threshold exceeded (",round(svalu[j],3),">",round(upper,3),").")
            }     
          }
        }
        if (sides==2) {
          if ((svalu[j]>=lower&svalu[j]<=upper)|is.na(svalu[j])) {
            valu[which(comp==compound)][j]<-svalu[j]
            note[which(comp==compound)][j]<-NA
          } else {
            valu[which(comp==compound)][j]<-NA
            if (class(x)!="data.frame") {
              valu[j]<-NA
            }
            if (svalu[j]<lower&!is.na(svalu[j])) {
              note[which(comp==compound)][j]<-paste0("Outlier: lower threshold not achieved (",round(svalu[j],3),"<",round(lower,3),").")
            }
            if (svalu[j]>upper&!is.na(svalu[j])) {
              note[which(comp==compound)][j]<-paste0("Outlier: upper threshold exceeded (",round(svalu[j],3),">",round(upper,3),").")
            }     
          }
        }
      } else {
        valu[which(comp==compound)][j]<-svalu[j]
        note[which(comp==compound)][j]<-NA
      }
      
    }
    glower[which(compos==compound)]<-lower
    gupper[which(compos==compound)]<-upper
    if (class(x)!="data.frame") {
      glower<-lower
      gupper<-upper
    }
    
    
    # Plotting.
    if (plot==TRUE) {
      
      # Stops in vector mode, when y is not compatible with x.
      if (class(x)!="data.frame") {
        if (length(x)!=length(y)) {
          stop(paste0("The length of vectors of concentrations and dates differ."))
        }
      }
      
      # If date is given only as the year, shifts the value to 1st june.
      if (max(nchar(as.character(sdate_start)),na.rm=TRUE)==4) {
        sdate_start<-paste0(x,"-06-01")
      }
      if (max(nchar(as.character(sdate_end)),na.rm=TRUE)==4) {
        sdate_end<-paste0(x,"-06-01")
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
      
      # Combines date_start and date_end to one date.
      sdate_start<-as.Date(as.character(sdate_start))
      sdate_end  <-as.Date(as.character(sdate_end))
      sdate      <-as.Date((as.numeric(sdate_start)+as.numeric(sdate_end))/2,origin="1970-01-01")
      
      valid<-which(!is.na(svalu)&!is.na(sdate))
      svalu<-svalu[valid]
      sdate<-sdate[valid]
      
      plot(as.Date(sdate,origin="1970-01-01"),svalu,cex=0,xlab=xlab,ylab=ylab,main=plotmain,pch=pch,ylim=c(1.1*min(svalu)-0.1*max(svalu),1.1*max(svalu)-0.1*min(svalu)))
      points(as.Date(sdate,origin="1970-01-01"),svalu,col=col.points,pch=pch)
      
      # Plotting oultiers
      if (length(svalu)>1) {
        polygon(c(min(sdate)-0.2*(max(sdate)-min(sdate)),min(sdate)-0.2*(max(sdate)-min(sdate)),max(sdate)+0.2*(max(sdate)-min(sdate)),max(sdate)+0.2*(max(sdate)-min(sdate))),c(1.2*max(svalu)-0.2*min(svalu),upper,upper,1.2*max(svalu)-0.2*min(svalu)),border=NA,col=rgb(0.6,0.3,0.3,0.3))
        points(as.Date(sdate[which(svalu>upper)],origin="1970-01-01"),svalu[which(svalu>upper)],cex=1.2,col="red",pch=19)
        if (sides==2) {
          polygon(c(min(sdate)-0.2*(max(sdate)-min(sdate)),min(sdate)-0.2*(max(sdate)-min(sdate)),max(sdate)+0.2*(max(sdate)-min(sdate)),max(sdate)+0.2*(max(sdate)-min(sdate))),c(1.2*min(svalu)-0.2*max(svalu),lower,lower,1.2*min(svalu)-0.2*max(svalu)),border=NA,col=rgb(0.3,0.3,0.6,0.3))
          points(as.Date(sdate[which(svalu<lower)],origin="1970-01-01"),svalu[which(svalu<lower)],cex=1.2,col="blue",pch=19)
        }
      }
    }
  }
  
  ## Generation of desired data type.
  # Output in genasis mode.
  if (!is.na(output)&output=="genasis") {
    result<-data.frame(valu,comp,date_start,date_end,temp,wind,note)
  }
  
  # Output in openair mode.
  if (!is.na(output)&output=="openair") {
    unidates<-as.data.frame(matrix(NA,0,2))
    for (a in unique(date_start)) {
      for (b in unique(date_end[which(date_start==a)])) {
        unidates<-rbind(unidates,as.Date(c(a,b),origin="1970-01-01"))
      }
    }
    
    unidates[,1]<-as.Date(unidates[,1],origin="1970-01-01")
    unidates[,2]<-as.Date(unidates[,2],origin="1970-01-01")
    
    result<-as.data.frame(matrix(NA,nrow(unidates),length(compos)+5))
    colnames(result)<-c("date","date_end","temp","wind","note",compos)
    

    for (i in 1:nrow(unidates)) {
      result[i,1]<-unidates[i,1]
      result[i,2]<-unidates[i,2]
      result[i,3]<-temp[which(date_start==unidates[i,1]&date_end==unidates[i,2]&comp==compound)][1]
      result[i,4]<-wind[which(date_start==unidates[i,1]&date_end==unidates[i,2]&comp==compound)][1]     
      if (paste(comp[which(date_start==unidates[i,1]&date_end==unidates[i,2]&!is.na(note))], collapse=', ')=="") {
        result[i,5]<-NA
      } else {
        result[i,5]<-paste0("Outliers: ",paste(comp[which(date_start==unidates[i,1]&date_end==unidates[i,2]&!is.na(note))], collapse=', '),".")
      }
      for (compound in compos) {
        result[i,compound]<-valu[which(date_start==unidates[i,1]&date_end==unidates[i,2]&comp==compound)][1]
      }
      if (length(which(date_start==unidates[i,1]&date_end==unidates[i,2]&comp==compound))>1) {
        warning(paste0("There was more than 1 record with the same start and end date for ",compound,", only the first record was used."))
        
      }
    }
    result[,1]<-as.Date(result[,1],origin="1970-01-01")
    result[,2]<-as.Date(result[,2],origin="1970-01-01")
  }
  
  # Output in vector-mode.
  if (class(x)!="data.frame"&is.na(output)) {
   result<-valu 
  }
  par(mfrow=c(1,1))
  return(list(res=result,lower=glower,upper=gupper))
}