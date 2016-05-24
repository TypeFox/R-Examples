formulaDeriv <-
function(inFormula, derivVar, data, nDerivs=1) {
  baseForm<-as.formula(inFormula)
  modMat <- lm(baseForm, data=data, x=T)$x
  baseLabels <- colnames(modMat)
  modLabels<-gsub("\\:","\\*",baseLabels,perl=T)
  modLabels <- sapply(modLabels, removeI)
  
  cleanLabels<-parse(text=modLabels)
  derivs<-vector("list",length(cleanLabels))
  nObs=dim(model.frame(baseForm,data=data))[1]
  outMat<-matrix(0, nrow=nObs, ncol=length(cleanLabels))
  
  if (baseLabels[1]=="(Intercept)") {
    for (i in 2:length(cleanLabels)) {     
      outMat[,i] <- matrix(0, nrow=nObs, ncol=1)
      derivs[[i]]<-cleanLabels[i]
      for (j in 1:nDerivs) {
        derivs[[i]]<-D(derivs[[i]],derivVar)
      }
      
      outMat[,i]<-matrix(eval(derivs[[i]], envir=data),nrow=nObs,ncol=1)
      
    }
  } else {
    for (i in 1:length(cleanLabels)) {
      
      derivs[[i]]<-cleanLabels[i]
      for (j in 1:nDerivs) {
        derivs[[i]]<-D(derivs[[i]], derivVar)
      }
      outMat[,i]<-matrix(eval(derivs[[i]], envir=data),nrow=nObs,ncol=1)
      
    }
  }
  return(outMat)
}
