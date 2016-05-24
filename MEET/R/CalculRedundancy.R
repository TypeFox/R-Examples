CalculRedundancy <-function(Factortrans,q,iicc){
	Prob<-as.numeric(iicc$background)
	Herror<-iicc$Herror
	HXmax<-iicc$HXmax
	
	probability_matrix<-probability(Factortrans,Prob)

	H<-switch(iicc$classentropy, "Shannon"=entropy.Shannon(probability_matrix),"Renyi"=entropy.Renyi(probability_matrix,q))

	Redundancia_corregida<-redundancy(H,HXmax,Herror,nrow(Factortrans),TRUE)
	
	return(Redundancia_corregida)
	
}

