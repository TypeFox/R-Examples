ModelEntropy<-function(iicc){
  require("MEET")	
  Prob<-iicc$background
  training.set<-Factortrans<-matriu<-iicc$Transcriptionfactor
    #iicc$missing.fun=TRUE
  q<-iicc$q<-iicc$correction_1rOrdre<-iicc$Herror<-iicc$HXmax<-iicc$Redundancia_corregida<-iicc$ErrorHX<-iicc$classentropy<-iicc$Entropy<-NULL
  q<-iicc$q<-iicc$parametersIdeal
 
  if (q==1) {iicc$classentropy<-"Shannon"
  }else{
			iicc$classentropy<-"Renyi"}
	
  iicc$correction_1rOrdre <- correction.entropy(q,nrow(matriu),1,iicc)
  iicc$Herror<-slot(iicc$correction_1rOrdre,"Herror")
  iicc$HXmax<-iicc$Herror[nrow(Factortrans)] 	
  iicc$Redundancia_corregida<-CalculRedundancy(matriu,q,iicc)
  iicc$ErrorHX<- iicc$ErrorHX<-slot(iicc$correction_1rOrdre,"sderror")[nrow(Factortrans)]

  memory<-Hmemory(iicc,Factortrans)
  iicc<-c(iicc,memory)
  parametersModel<-list(Order=iicc$q,HXmax=iicc$HXmax,ErrorHX=iicc$ErrorHX,Redundancy=iicc$Redundancia_corregida)		
  y<-list(model=iicc$Entropy,parameterModel=parametersModel)	
  return(y)	


}
