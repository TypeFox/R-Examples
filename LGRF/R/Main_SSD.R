
# Some functions used in this file are from the SKAT package

LGRF.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, result.null, Gsub.id=NULL, interGXT=FALSE, similarity='GR', impute.method='fixed', MinP.compare=FALSE, ...){
  
  id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
  if(length(id1) == 0){
    MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
    stop(MSG)
  }	
  SetID<-SSD.INFO$SetInfo$SetID[id1]
  
  try1<-try(Get_Genotypes_SSD(SSD.INFO, SetIndex, is_ID=T),silent = TRUE)
  if(class(try1) != "try-error"){
    Z<-try1
    Is.Error<-FALSE	
  } else {
    err.msg<-geterrmessage()
    msg<-sprintf("Error to get genotypes of %s: %s",SetID, err.msg)
    stop(msg)
  }
  #re1 = LGRF.SSD.GetSNP_Weight(SSD.INFO, SetIndex, obj.SNPWeight=NULL)
  re<-test.LGRF(Z, Gsub.id=Gsub.id, result.null, interGXT=interGXT,similarity=similarity,impute.method=impute.method)
  # compare with MinP test
  if(MinP.compare==TRUE){
    re.MinP<-test.MinP(Z, Gsub.id=Gsub.id, result.null, impute.method=impute.method, ...)
    re$p.value<-cbind(re$p.value,re.MinP$p.value);colnames(re$p.value)<-c('LGRF','MinP')
  }
  
  return(re)
}

#LGRF.SSD.All function
LGRF.SSD.All = function(SSD.INFO, result.null, Gsub.id=NULL, interGXT=FALSE, similarity='GR', impute.method='fixed', MinP.compare=FALSE, ...){
  N.Set<-SSD.INFO$nSets
  OUT.Pvalue<-rep(NA,N.Set)
  OUT.Pvalue.MinP<-rep(NA,N.Set)
  OUT.Marker<-rep(NA,N.Set)
  #OUT.Marker.Test<-rep(NA,N.Set)
  OUT.Error<-rep(-1,N.Set)
  
  for(i in 1:N.Set){
    Is.Error<-TRUE
    try1 = try(LGRF.SSD.OneSet_SetIndex(SSD.INFO=SSD.INFO, SetIndex=i, result.null=result.null, Gsub.id=Gsub.id, interGXT=interGXT, MinP.compare=MinP.compare, ...))
    
    if(class(try1) != "try-error"){
      re<-try1; print(SSD.INFO$SetInfo[i,1:3])
      Is.Error<-FALSE
    } else {
      
      err.msg<-geterrmessage()
      msg<-sprintf("Error to run LGRF for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
      warning(msg,call.=FALSE)
      
    }
    
    if(!Is.Error){
      
      OUT.Pvalue[i]<-re$p.value[1]
      if(MinP.compare==TRUE){OUT.Pvalue.MinP[i]<-re$p.value[2]}
      OUT.Marker[i]<-re$n.marker
      #OUT.Marker.Test[i]<-re$param$n.marker.test
    }
  }
  
  out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue, N.Marker=OUT.Marker)
  if(MinP.compare==TRUE){out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue,P.value_MinP=OUT.Pvalue.MinP, N.Marker=OUT.Marker)}
  re<-list(results=out.tbl)
  class(re)<-"LGRF_SSD_ALL"
  
  return(re)
}

