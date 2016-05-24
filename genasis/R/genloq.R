genloq<-function(x,y=NA,input="openair",output=NA,loq=NA,method="mle",distr="lnorm",pollutant=NA,plot=TRUE,ylim=c(NA,NA),columns=2,col.points="black",pch=1,xlab="Date",ylab="Concentration",main=NA) {
  
  ## "Genmle" function definition.
  genmle<-function(tvalu,distr,loq) {
    
    if (length(which(is.na(suppressWarnings(as.numeric(tvalu)))))>0) {
      
      # Prepares numeric vector "nvalu" and estimates loq.
      suppressWarnings(nvalu<-as.numeric(gsub(",",".",tvalu)))
      if (distr=="lnorm") {
        nvalu<-log(nvalu)
        loq<-log(loq)
      }
      
      # Creates a data frame "field" with censoring intervals.
      nvalusloq<-nvalu
      nvalusloq[which(is.na(nvalusloq))]<-loq
      field<-data.frame(nvalu,nvalusloq)
      colnames(field)<-c("left","right")
      
      # Finds a best-fitting normal distribution.
      nloq<-length(which(tvalu=="loq"|tvalu=="LoQ"|tvalu=="LOQ"))
      mmm<-fitdistcens(field,distr="norm")$estimate[1]
      sss<-fitdistcens(field,distr="norm")$estimate[2]
      
      # Define by-magnitude ordered values for loq substitution.
      substituents<-qnorm((1:nloq)/(nloq+1)*pnorm(loq,mmm,sss),mmm,sss)
      
      # Computes order of loq values and substitute them.
      qvalu<-nvalu
      distance<-c(NA,rep(Inf,length(tvalu)),NA)
      distance[which(!is.na(tvalu)&tvalu!="loq"&tvalu!="LOQ"&tvalu!="LoQ")+1]<-0
      for (i in 1:(length(distance)/2)) {
        for (j in 2:(length(distance)-1)) {
          distance[j]<-min(distance[j-1]+1,distance[j],distance[j+1]+1,na.rm=TRUE)
        }
      }
      distance<-distance[c(-1,-length(distance))]
      
      for (j in 1:max(distance,na.rm=TRUE)) {
        thesearej<-which(distance==j)
        averagesj<-rep(NA,length(thesearej))
        for (k in 1:length(averagesj)) {
          averagesj[k]<-mean(c(qvalu[thesearej[k]-1],qvalu[thesearej[k]+1]),na.rm=TRUE)
        }
        qvalu[thesearej[order(averagesj)]]<-substituents[(length(substituents)-length(thesearej)+1):length(substituents)]
        substituents<-substituents[1:(length(substituents)-length(thesearej))]
      }
      
      # Combines original "valu" and quasiresult "qvalu" to a final result "rvalu".
      if (distr=="lnorm") {
        qvalu<-exp(qvalu)
      }    
    } else {
      qvalu<-as.numeric(tvalu)
    }   
    return(qvalu)
  }
  # End of "genmle" function definition.
  
  
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
    valu      <-gsub(",",".",as.character(x[,1]))
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
      valu      <-c(valu,as.character(gsub(",",".",(x[,compound]))))
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
    valu      <-as.character(gsub(",",".",(x)))
    comp      <-rep(pollutant,length(x))
    date_start<-y
    date_end  <-y
    temp      <-rep(NA,length(x))
    wind      <-rep(NA,length(x)) 
  }
  
  # Settings of plotting output.
  if (plot==TRUE&length(compos)>1) {
    par(mfrow=c(ceiling(length(compos)/columns),columns))
  }
  if (plot==TRUE&length(compos)<=1) {
    par(mfrow=c(1,1))
  }
  if (plot==TRUE&length(compos)>6) {
    warning(paste0("There is a lot of compounds to plot in the dataset x (",length(compos),"). The plot could be crowded on some devices."))
  }
  
  ## Verification and sorting.
  valid<-which(!is.na(date_start)|!is.na(date_end))
  
  valu      <-valu[valid]
  comp      <-comp[valid]
  date_start<-date_start[valid]
  date_end  <-date_end[valid]
  temp      <-temp[valid]
  wind      <-wind[valid]
  
  date<-as.Date((as.numeric(gendate(date_start))+as.numeric(gendate(date_end)))/2,origin="1970-01-01")
  
  valu      <-valu[order(date)]
  comp      <-comp[order(date)]
  date_start<-date_start[order(date)]
  date_end  <-date_end[order(date)]
  temp      <-temp[order(date)]
  wind      <-wind[order(date)]
  date      <-date[order(date)]
  
  
  ylim<-as.numeric(ylim)
  
  ## Loq substitution
  note<-rep(NA,length(valu))
  gloq<-c()
  
  for (compound in compos) {
    svalu<-valu[which(comp==compound)]
    sdate_start<-date_start[which(comp==compound)]
    sdate_end<-date_end[which(comp==compound)]
    sdate<-date[which(comp==compound)]
    
    if (class(x)!="data.frame") {
      svalu<-x[order(date)]
      sdate_start<-y[order(date)]
      sdate_end<-y[order(date)]
      sdate<-date[order(date)]
    }
    
    # Find "loq" if not set explicitely.
    if (is.na(loq)) {
      sloq<-suppressWarnings(min(as.numeric(gsub(",",".",svalu)),na.rm=TRUE))
    } else {
      sloq<-loq
    }
    
    # Saves positions of NAs in valu and creates "ssvalu".
    notna<-which(!is.na(svalu))
    ssvalu<-svalu[notna]
    
    qvalu<-suppressWarnings(as.numeric(gsub(",",".",ssvalu)));
    toreplace<-which(is.na(qvalu))
    
    if (method=="mle") {qvalu<-genmle(ssvalu,distr=distr,loq=sloq)}
    if (method=="exc") {}
    if (method=="2.0") {qvalu[toreplace]<-sloq/2}
    if (method=="1.0") {qvalu[toreplace]<-sloq}
    if (method=="sq2") {qvalu[toreplace]<-sloq/sqrt(2)}
    
    rvalu<-suppressWarnings(as.numeric(gsub(",",".",svalu)))
    rvalu[notna]<-qvalu
    
    wasloq<-notna[toreplace]
    
    valu[which(comp==compound)]<-rvalu
    note[which(comp==compound)][wasloq]<-"LoQ treated."
    
    # Plotting.
    if (plot==TRUE) {
      
      # Stops in vector mode, when y is not compatible with x.
      if (class(x)!="data.frame") {
        if (length(x)!=length(y)) {
          stop(paste0("The length of vectors of concentrations and dates differ."))
        }
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
      
      # Plot.
      if (!is.na(max(ylim))) {
        plotylim<-ylim
      } else {
        plotylim<-c(min(rvalu,na.rm=TRUE),max(rvalu,na.rm=TRUE))
      }
      
      plot(sdate,rvalu,cex=0,ylim=plotylim,xlab=xlab,ylab=ylab,main=plotmain,pch=pch)
      points(sdate,suppressWarnings(as.numeric(gsub(",",".",svalu))),col=col.points,pch=pch)
      
      # Plotting LoQs and their area.
      polygon(c(min(sdate)-0.2*(max(sdate)-min(sdate)),min(sdate)-0.2*(max(sdate)-min(sdate)),max(sdate)+0.2*(max(sdate)-min(sdate)),max(sdate)+0.2*(max(sdate)-min(sdate))),c(1.2*min(rvalu,na.rm=TRUE)-0.2*max(rvalu,na.rm=TRUE),sloq,sloq,1.2*min(rvalu,na.rm=TRUE)-0.2*max(rvalu,na.rm=TRUE)),border=NA,col=rgb(0.3,0.3,0.6,0.3))
      points(sdate[wasloq],rvalu[wasloq],cex=1.2,col="blue",pch=19)
      }
    
    gloq<-c(gloq,sloq)
    }
  
  valu<-as.numeric(valu)

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
        result[i,5]<-paste0("LoQs treated for: ",paste(comp[which(date_start==unidates[i,1]&date_end==unidates[i,2]&!is.na(note))], collapse=', '),".")
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
    result<-rvalu 
  }
  par(mfrow=c(1,1))
  return(list(res=result,loq=gloq))
}