methStatusEval <- function(x, error=0.05, uninformative=TRUE){
	#Evaluate for classification as 'methylation-susceptible locus' or 'non methylated"
	#Uses mixed data: 11, 10, 1, 0 are the four possible states for both enzymes
	threshold= ((error + error) - (2*error*error))
	#A locus is classified as "methylation-susceptible" if the frecuency of methylation evidence (10, 1 or 0) is above threshold
	#If 00 pattern assumes uninformative state. 
	if(uninformative) fMet <- 1-(length(x[which(x==11|x==0)])/length(x)) 
	# If 00 patter assumes hypermethylation 
	else fMet <- 1-(length(x[which(x==11)])/length(x))
	return(fMet>threshold)
}
