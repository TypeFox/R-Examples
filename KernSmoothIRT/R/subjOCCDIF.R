subjOCCDIF <- function(x, stype = c("ObsScore","Theta","ThetaML","ScoreML")){

	if(missing(stype)) stype <- "ObsScore"
	OCCs <- list(length(x$DIF))
	
	for (i in 1:length(x$DIF)){
		OCCs[[i]] <- subjOCC(x$DIF[[i]],stype=stype)
	}
	
	return(OCCs)
	
}