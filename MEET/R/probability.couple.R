probability.couple <-function(Prob){
	out<-.C("probabilitycouple",background=as.double(Prob),Probtrans=double(16))
	return(out$Probtrans)
	}

