

	getLik<-function(x,answers){
		

		check<-numeric()

		for(i in 1:length(x)){
			if(answers[i]==0){check[i]<-(1-x[i])}
			else{check[i]<-x[i]}
		}			

		return(prod(check))
	

	}


