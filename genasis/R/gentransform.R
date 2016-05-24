gentransform<-function(x,y=NA,input="openair",output=NA,pollutant=NA) {
  
  ## Selection of relevant compound(s).
  # If x is defined as genasis-type data frame, load the compounds.
  if (class(x)=="data.frame"&input=="genasis") {
    compounds<-unique(x[,2])
  }
  
  # If x is defined as openair-type data frame, load the compounds.
  if (class(x)=="data.frame"&input=="openair") {
    compounds<-colnames(x)[-which(is.element(colnames(x),c("date","date_end","temp","wind")))]
  }
  
  # If x is defined as numeric vector, set compounds=pollutant.
  if (class(x)!="data.frame") {
    compounds<-pollutant
  }
  
  compounds<-as.character(compounds)
  
  
  ## Creating column variables.
  
  # If x is defined as genasis-type data frame, set column vectors.
  if (class(x)=="data.frame"&input=="genasis") {
    valu      <-as.numeric(x[,1])
    comp      <-as.character(x[,2])
    date_start<-as.Date(as.character(x[,3]))
    date_end  <-as.Date(as.character(x[,4]))
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
    for (compound in compounds) {
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
  
  ## Generation of desired data type.
  if (output=="genasis") {
    result<-data.frame(valu,comp,date_start,date_end,temp,wind)
  }
  
  if (output=="openair") {
    unidates<-as.data.frame(matrix(NA,0,2))
    for (a in unique(date_start)) {
      for (b in unique(date_end[which(date_start==a)])) {
        unidates<-rbind(unidates,as.Date(c(a,b),origin="1970-01-01"))
      }
    }
    
    unidates[,1]<-as.Date(unidates[,1],origin="1970-01-01")
    unidates[,2]<-as.Date(unidates[,2],origin="1970-01-01")
    
    result<-as.data.frame(matrix(NA,nrow(unidates),length(compounds)+4))
    colnames(result)<-c("date","date_end","temp","wind",compounds)
    

    for (i in 1:nrow(unidates)) {
      result[i,1]<-unidates[i,1]
      result[i,2]<-unidates[i,2]
      result[i,3]<-temp[which(date_start==unidates[i,1]&date_end==unidates[i,2])][1]
      result[i,4]<-wind[which(date_start==unidates[i,1]&date_end==unidates[i,2])][1]
      for (compound in compounds) {
        result[i,compound]<-valu[which(date_start==unidates[i,1]&date_end==unidates[i,2]&comp==compound)][1]
      }
      if (length(which(date_start==unidates[i,1]&date_end==unidates[i,2]&comp==compound))>1) {
        warning(paste0("There was more than 1 record with tha same start and end date for ",compound,", only the first record was used."))
        
      }
    }
    result[,1]<-as.Date(result[,1],origin="1970-01-01")
    result[,2]<-as.Date(result[,2],origin="1970-01-01")
  }
  return(res=result)
}