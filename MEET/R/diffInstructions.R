diffInstructions <-function (training.set,HX,HXmax,Herror,Redundancia_corregida){	
	
		Rsequencia<-redundancy(HX,HXmax,Herror,finite=nrow(training.set),TRUE)
		x<-sum(Redundancia_corregida*(abs(Rsequencia-Redundancia_corregida)))
		    
		return(x)
}

