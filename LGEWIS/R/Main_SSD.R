
# Some functions used in this file are from the SKAT package

## Marginal genetic association

GA.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, result.prelim, Gsub.id=NULL, MinP.adjust=NULL, ...){

  id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
  if(length(id1) == 0){
    MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
    stop(MSG)
  }
  SetID<-SSD.INFO$SetInfo$SetID[id1]

  try1<-try(Get_Genotypes_SSD(SSD.INFO, SetIndex, is_ID=T),silent = TRUE)
  if(class(try1) != "try-error"){
    G<-try1
    Is.Error<-FALSE
  } else {
    err.msg<-geterrmessage()
    msg<-sprintf("Error to get genotypes of %s: %s",SetID, err.msg)
    stop(msg)
  }
  re<-GA.test(result.prelim, G, Gsub.id=Gsub.id, MinP.adjust=MinP.adjust, ...)

  return(re)
}

#LGRF.SSD.All function
GA.SSD.All = function(SSD.INFO, result.prelim, Gsub.id=NULL, MinP.adjust=NULL, ...){
  N.Set<-SSD.INFO$nSets
  OUT.Pvalue<-rep(NA,N.Set)
  OUT.Pvalue.MinP<-rep(NA,N.Set)
  OUT.Marker<-rep(NA,N.Set)
  #OUT.Marker.Test<-rep(NA,N.Set)
  OUT.Error<-rep(-1,N.Set)
  OUT.single<-c()

  for(i in 1:N.Set){
    if(i%%100==0){print(paste0(i," sets finished"))}
    Is.Error<-TRUE
    try1 = try(GA.SSD.OneSet_SetIndex(SSD.INFO=SSD.INFO, SetIndex=i, result.prelim=result.prelim, Gsub.id=Gsub.id, MinP.adjust=MinP.adjust, ...))

    if(class(try1) != "try-error"){
      re<-try1;
      Is.Error<-FALSE
    } else {

      err.msg<-geterrmessage()
      msg<-sprintf("Error to run GA for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
      warning(msg,call.=FALSE)

    }
    print(paste(i,"-th set",SSD.INFO$SetInfo[i,2],",",SSD.INFO$SetInfo[i,3],"SNPs"))

    if(!Is.Error){

      OUT.Pvalue[i]<-re$p.value
      OUT.Pvalue.MinP[i]<-re$p.MinP
      OUT.Marker[i]<-re$n.marker
      temp.single<-cbind(rep(SSD.INFO$SetInfo[i,2],nrow(re$p.single)),rownames(re$p.single),re$p.single)
      rownames(temp.single)<-NULL
      OUT.single<-rbind(OUT.single,temp.single)

      #OUT.Marker.Test[i]<-re$param$n.marker.test
    }
  }

  out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue, N.Marker=OUT.Marker)
  out.tbl.single<-data.frame(OUT.single);colnames(out.tbl.single)<-c('Region.name','SNP.name','MAF','p.value')

  if(length(MinP.adjust)!=0){out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue,P.value_MinP=OUT.Pvalue.MinP, N.Marker=OUT.Marker)}
  re<-list(results.single=out.tbl.single,results=out.tbl)
  class(re)<-"GA_SSD_ALL"

  return(re)
}

## Gene-environement intereaction

GEI.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, result.prelim, Gsub.id=NULL, MinP.adjust=NULL, ...){

  id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
  if(length(id1) == 0){
    MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
    stop(MSG)
  }
  SetID<-SSD.INFO$SetInfo$SetID[id1]

  try1<-try(Get_Genotypes_SSD(SSD.INFO, SetIndex, is_ID=T),silent = TRUE)
  if(class(try1) != "try-error"){
    G<-try1
    Is.Error<-FALSE
  } else {
    err.msg<-geterrmessage()
    msg<-sprintf("Error to get genotypes of %s: %s",SetID, err.msg)
    stop(msg)
  }
  re<-GEI.test(result.prelim, G, Gsub.id=Gsub.id, MinP.adjust=MinP.adjust, ...)

  return(re)
}

#GEI.SSD.All function
GEI.SSD.All = function(SSD.INFO, result.prelim, Gsub.id=NULL, MinP.adjust=NULL, ...){
  N.Set<-SSD.INFO$nSets
  OUT.Pvalue<-rep(NA,N.Set)
  OUT.Pvalue.MinP<-rep(NA,N.Set)
  OUT.Marker<-rep(NA,N.Set)
  #OUT.Marker.Test<-rep(NA,N.Set)
  OUT.Error<-rep(-1,N.Set)
  OUT.single<-c()

  for(i in 1:N.Set){
    if(i%%100==0){print(paste0(i," sets finished"))}
    Is.Error<-TRUE
    #try1 = try(GEI.SSD.OneSet_SetIndex(SSD.INFO=SSD.INFO, SetIndex=i, result.prelim=result.prelim, Gsub.id=Gsub.id, MinP.adjust=MinP.adjust, ...))
    try1 = try(GEI.SSD.OneSet_SetIndex(SSD.INFO=SSD.INFO, SetIndex=i, result.prelim=result.prelim, Gsub.id=Gsub.id, MinP.adjust=MinP.adjust, ...))

    if(class(try1) != "try-error"){
      re<-try1;
      Is.Error<-FALSE
    } else {

      err.msg<-geterrmessage()
      msg<-sprintf("Error to run GEI for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
      warning(msg,call.=FALSE)

    }
    print(paste(i,"-th set",SSD.INFO$SetInfo[i,2],",",SSD.INFO$SetInfo[i,3],"SNPs"))

    if(!Is.Error){

      OUT.Pvalue[i]<-re$p.value
      OUT.Pvalue.MinP[i]<-re$p.MinP
      OUT.Marker[i]<-re$n.marker
      temp.single<-cbind(rep(SSD.INFO$SetInfo[i,2],nrow(re$p.single)),rownames(re$p.single),re$p.single)
      rownames(temp.single)<-NULL
      OUT.single<-rbind(OUT.single,temp.single)

      #OUT.Marker.Test[i]<-re$param$n.marker.test
    }
  }

  out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue, N.Marker=OUT.Marker)
  out.tbl.single<-data.frame(OUT.single);colnames(out.tbl.single)<-c('Region.name','SNP.name','MAF','p.value')

  if(length(MinP.adjust)!=0){out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue,P.value_MinP=OUT.Pvalue.MinP, N.Marker=OUT.Marker)}
  re<-list(results.single=out.tbl.single,results=out.tbl)
  class(re)<-"GEI_SSD_ALL"

  return(re)
}



