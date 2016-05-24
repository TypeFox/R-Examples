chooseModel<-function(AUC, iicc){
 require("MEET")

    variance<-unlist(lapply(AUC, function(x){sd(x)}))
    mean<-unlist(lapply(AUC, function(x){mean(x)}))
    areamv<-mean*(1-variance)
    indexParameter<-which(areamv==max(areamv))
    parameterIdeal<-iicc$vector[indexParameter]
	y<-list(parameterIdeal=parameterIdeal,index=indexParameter)
    return(y)
     
}

