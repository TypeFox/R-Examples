ABdescrip <-
  function(behavior,PhaseX){
    graphics.off()
    boxplot(behavior~PhaseX)
    ab<-NULL
    ab<<-recordPlot()
    writeLines("-----------n-------------")
    t1<-table(PhaseX)
    print(t1)
    abmean<-tapply(behavior,PhaseX,mean,na.rm=T)
    pmean<-c(round(abmean,3))
    writeLines("-----------mean-------------")
    print(pmean)
    
    tabmean<-tapply(behavior,PhaseX,mean,trim=.1,na.rm=T)
    tmean<-c(round(tabmean,3))
    writeLines("-----------10% trim mean-------------")
    print(tmean)
    
    abmedian<-tapply(behavior,PhaseX,median,na.rm=T)
    pmedian<-c(round(abmedian,3))
    writeLines("----------median------------")
    print(pmedian)
    absd<-tapply(behavior,PhaseX,sd,na.rm=T)
    psd<-c(round(absd,3))
    writeLines("------------SD--------------")
    print(psd)
    cv<-psd/pmean
    pcv<-c(round(cv,3))
    writeLines("------------CV--------------")
    print(pcv)
    
    writeLines("---------range----------")
    abrange<-tapply(behavior,PhaseX,range,na.rm=T)
    prange<- do.call("rbind", abrange)
    print(prange)
    
    
    writeLines("---------iqr----------")
    abiqr<-tapply(behavior,PhaseX,IQR,na.rm=T)
    print(abiqr)
    
    writeLines("---------quantiles----------")
    abquan<-tapply(behavior,PhaseX,quantile,na.rm=T)
    do.call("rbind", abquan)
    
    
  }
