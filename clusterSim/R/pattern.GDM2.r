.manualP<-function(patternManual,data,nomOptValues){
  ret<-NULL
  for(i in 1:length(patternManual)){
    p<-patternManual[i]
    if(p=="min") {
      ret<-c(ret,min(data[,i]))
    }
    else if(p=="max") {
      ret<-c(ret,max(data[,i]))
    }
    else if(p=="nom"){
      ret<-c(ret,nomOptValues[i])
    }
    else{
      ret<-c(ret,as.numeric(p))
    }
    }
  ret
}

pattern.GDM2<-function(data,performanceVariable,nomOptValues=NULL,weightsType="equal",weights=NULL,patternType="upper",patternCoordinates="dataBounds",patternManual=NULL,nominalTransfMethod=NULL){
#WCZYTANIE DANYCH
pdata<-data
data<-as.matrix(data)

#DEFINICJE PARAMETROW

#SPRAWDZENIE PARAMETROW
for(v in performanceVariable){
  if(sum(v==c("s","n","d"))==0) stop("performance variable should be one of the following: s-stimuli,d-destimuli or n-nominant variable")
}
if (length(performanceVariable)!=ncol(data)) stop("performance variable vector should have the size equal to numberof variables")
if(sum(performanceVariable=="n")==0)nomOptValues<-rep(0,ncol(data))
if (length(nomOptValues)!=ncol(data)) stop("vector of optimal values for nominant variables should have the size equal to numberof variables")
#if (length(noCategories)!=ncol(data)) stop("Number of categories vector should have the size equal to numberof variables")
if(sum(patternType==c("upper","lower"))==0) stop ("pattern should be one of the following:  \"upper\",\"lower\"")
if(sum(patternCoordinates==c("dataBounds","bounds","manual"))==0) stop ("pattern should be one of the following:  \"dataBounds\",\"bounds\",\"manual\"")
if (patternCoordinates=="manual" && length(patternManual)!=ncol(data)) stop("pattern manual value vector should have the size equal to number of variables")
if(sum(performanceVariable=="n")==0){
#print("zmienilem")
nominalTransfMethod<-rep("symmetrical",ncol(data))
}
if(is.null(nominalTransfMethod))nominalTransfMethod=rep("none",ncol(data))
if(length(nominalTransfMethod)==1)nominalTransfMethod=rep(nominalTransfMethod,ncol(data))
if(length(nominalTransfMethod)!=ncol(data))stop("vector of transfer methods for nominant variables should have the size equal to numberof variables")
for(v in nominalTransfMethod){
  if(sum(v==c("symmetrical","database","none"))==0) stop("transfer methods should be one of the following: s-symmetrical,d-database")
}

#ZNALEZIENIE WZORCA
pattern<-rep(0,ncol(data))
vTypes<-performanceVariable
for(i in 1:length(performanceVariable)){
  if(patternType=="lower"){
  if(performanceVariable[i]=="s"){
    vTypes[i]<-"d"
  }
  if(performanceVariable[i]=="d"){
    vTypes[i]<-"s"
  }
  if(performanceVariable[i]=="n"){
    vTypes[i]<-"a"
  }
  if(performanceVariable[i]=="a"){
    vTypes[i]<-"n"
  }
  }
}
for(i in 1:length(vTypes)){
      #if(patternCoordinates=="bounds" && vTypes[i]=="s"){
      #  pattern[i]<-noCategories[i]
      #}
      #if(patternCoordinates=="bounds" && vTypes[i]=="d"){
      #  pattern[i]<-1
      #}
      if(patternCoordinates=="dataBounds" && vTypes[i]=="s"){
        pattern[i]<-max(data[,i])
      }
      if(patternCoordinates=="dataBounds" && vTypes[i]=="d"){
        pattern[i]<-min(data[,i])
      }
      if(vTypes[i]=="n" ){
        pattern[i]<-nomOptValues[i]
      }
      if(vTypes[i]=="a"){
        if(nominalTransfMethod[i]!="none"){
        if(nominalTransfMethod[i]=="symmetrical" ){
          #print("symmetrical")
          t<-c(nomOptValues[i],data[,i])
          t1<-unique(t)
          dim(t1)<-c(length(t1),1)
          d<-as.matrix(GDM2(t1))[1,]
          for(j in 1:nrow(data)){
            for( k in 1:(length(d))){
            if(data[j,i]==t1[k]){
            data[j,i]<-d[k]
            }
            }
          }
          dd<-cbind(t1,d)
          nomOptValues[i]<-0
        }
        else{
          #print("database")
          t<-c(nomOptValues[i],data[,i])
          dim(t)<-c(nrow(data)+1,1)
          d<-as.matrix(GDM2(t))[1,]
          for(j in 1:nrow(data)){
            data[j,i]<-d[j+1]
          }
          dd<-d
          nomOptValues[i]<-0
       }
       }
     }
}

if(patternCoordinates=="manual"){
  pattern<-.manualP(patternManual,data,nomOptValues)
}

if(patternType=="lower"){
#print("Data after transformation")
#print(data)
}
#print(paste("Pattern ",paste(pattern,collapse=",")))

#WLASCIWE PORZADKOWANIE
cdata<-rbind(pattern,data)
gdm<-GDM2(cdata,weightsType=weightsType, weights=weights)
gdm_p<-as.matrix(gdm)[1,][-1]
names(gdm_p)<-row.names(data)
#print("GDM distances from pattern object")
#print(gdm_p)
#print("Sorted GDM distances form pattern object")
if(patternType=="upper"){
sortedgdm_p<-sort(gdm_p)
}
if(patternType=="lower"){
sortedgdm_p<-sort(gdm_p,decreasing=TRUE)
}
#print(sortedgdm_p)
resul<-list(pdata=pdata,data=rbind(data,pattern),distances=gdm_p,sortedDistances=sortedgdm_p)
resul
}
