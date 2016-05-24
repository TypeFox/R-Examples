posClusters <-
function(subTax,k){
	results<-list()
	cB<-comBim(subTax,k)
	for (i in 1:length(cB[,1])){
		VecBin<-cB[i,]
		results[[i]]<-SetUpGroupMatrix(subTax,VecBin)
		}
	res<-list()	
	for (i in 1:length(results)){
		res[[i]]<-SumSetUpGroupMatrix(results[[i]])
		}
	res<-RemRedun(res)	
	return(res)
	
	}

