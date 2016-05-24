ROCmodel <-function(Scores, labels, iicc){

	require("ROCR")
	
	AREA<-perf<-lapply(seq(1, length(Scores), 1), function(x){})
	 
	for (j in c(1:length(Scores))){
    
	  pred<-prediction(Scores[[j]],labels[[j]])
	  perf[[j]]<-performance(pred,"tpr","fpr")
	  r<-performance(pred,'auc')
	  AREA[[j]]<-as.numeric(slot(r,"y.values"))
	 }

  	bestparameter<-chooseModel(AREA, iicc)
	iicc$parameterIdeal<-iicc$parametersIdeal<-bestparameter$parameterIdeal
 	FinalModel<- Models(iicc)
 	iicc$model<-FinalModel
	y<-list(Area=AREA, parameters=iicc$parametersIdeal, model=iicc$model)
	
}

