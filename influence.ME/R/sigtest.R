sigtest <- function(estex, test=1.96, parameters=0, sort=FALSE, to.sort=NA) 
	{
	
		n.groups 					<- dim(estex$alt.test)[1]
		n.parameters 				<- dim(estex$alt.test)[2]
		ifelse(parameters==0, sel 	<- 1:n.parameters, sel <- parameters)	
		
		posneg 						<- ifelse(test>0, 1, -1)
		
		original.teststat 			<- as.matrix(estex$or.test[sel])
		original.test 				<- posneg * original.teststat > posneg * test
		
		altered.teststat 			<- as.matrix(estex$alt.test[, sel])
		altered.test				<- posneg * altered.teststat > posneg * test
		
		altered.change				<- t(apply(altered.test, 1, function(x) x!=original.test))
		colnames(altered.change)	<- colnames(altered.test)		
				
		output <- list()
		
		for(i in 1:n.parameters)
			{
			parname <- colnames(altered.test)[i]
			if(parname=="(Intercept)") parname<-"Intercept"
			
			output[[parname]] <- data.frame(Altered.Teststat=altered.teststat[,i], Altered.Sig=altered.test[,i], Changed.Sig=altered.change[,i])
			}
		
		if(sort == TRUE)
			{
			if(length(names(output)) == 1)
				{
				to.sort <- names(output)
				}
	
			if(length(names(output)) > 1 & is.na(to.sort))
				{
				stop("Indicate on which variable to sort the output: please specify to.sort")
				}
			
		
			ordering 	<- order(output[[to.sort]][,1])	
			output 		<- lapply(output, function(x) x[ordering,])
			}

	return(output)
                
	}


