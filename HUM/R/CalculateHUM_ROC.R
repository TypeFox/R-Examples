CalculateHUM_ROC<-function(data,indexF,indexClass,indexLabel,seq)
{
  
  indexL=NULL
  
  label=unique(data[,indexClass])
  for(i in 1:length(indexLabel))
  {
  indexL=c(indexL,which(label==indexLabel[i]))
  }
  
  Sp=NULL
  Sn=NULL
  S3=NULL
  optSp=NULL
  optSn=NULL
  optS3=NULL
  
  indexEach=NULL
  indexUnion=NULL
    
  for(i in 1:length(label))
  {
      vrem=which(data[,indexClass]==label[i])
      indexEach=c(indexEach,list(vrem))
      if(label[i]%in%indexLabel)
      indexUnion=union(indexUnion,vrem)
  }
  
  for(i in 1:length(indexF))
  {
    s_data=NULL
    dataV=data[,indexF[i]]	
    prodValue=1
    for (j in 1:length(indexLabel))
    {
    vrem=sort(dataV[indexEach[[indexL[j]]]])
        
    s_data=c(s_data,list(vrem))
    prodValue = prodValue*length(vrem)
    }
    #claculate the threshold values for plot of 2D ROC and 3D ROC
    thresholds <- sort(unique(dataV[indexUnion]))
    thresholds=(c(-Inf, thresholds) + c(thresholds, +Inf))/2
    
    out=CalcROC(s_data,seq[,indexF[i]], thresholds)

    Sp=c(Sp,list(out$Sp))
    Sn=c(Sn,list(out$Sn))
    
    optSp=c(optSp,out$optSp)
    optSn=c(optSn,out$optSn)
    
    if(length(indexLabel)==3)
    {
      S3=c(S3,list(out$S3))
      optS3=c(optS3,out$optS3)
    }
  
  }
  
  names(optSp)=names(data[,indexF,drop = FALSE])
  names(optSn)=names(data[,indexF,drop = FALSE])
  
  names(Sp)=names(data[,indexF,drop = FALSE])
  names(Sn)=names(data[,indexF,drop = FALSE])
  if(length(indexLabel)==3)
  {
  names(S3)=names(data[,indexF,drop = FALSE])
  names(optS3)=names(data[,indexF,drop = FALSE])
  }
  
  if(length(indexLabel)==3)
  {
  return(list(Sp=Sp,Sn=Sn,S3=S3,thresholds=thresholds,optSn=optSn,optSp=optSp,optS3=optS3))
  }
  else
  {
    return(list(Sp=Sp,Sn=Sn,optSn=optSn,optSp=optSp))
  }
}
