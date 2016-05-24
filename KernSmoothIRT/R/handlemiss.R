handlemiss <-
function(x,miss){
	## This function takes a single subjects' responses and treats missing values accordingly.
	answers<-x

	


	if(miss=="random.multinom"){
		## For each missing item, sample from a multinomial distribution with probabilities equal to the relative frequencies of each option amongst
		## non-missing responses.
		freqbycol<-table(answers)
		
		answers[is.na(answers)]<-as.numeric(names(freqbycol))[which(rmultinom(1,1,freqbycol)==1)]


	}
	
	else if(miss=="random.unif"){
		## For each missing item, sample from a discrete normal distribution with probabilities equal probabilities for all options.
		answers[is.na(answers)]<-sample(na.omit(unique(answers)),1)

	}
	

	else if(miss=="option"){
		## Treat missing values as categorical. Use -1 as an identifier instead of NA.
	
		answers[is.na(answers)]=-1

	
	}
	


	return(answers)



}

