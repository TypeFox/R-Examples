lrtCapwire <-
function(ecm, tirm, bootstraps=100){
	
	likE <- ecm$likelihood
	
	likT <- tirm$likelihood
	
	LR <- likT - likE # Test Statistic
	
	# Simulate data sets under ECM using MLE of N from ECM model
	
	n <- ecm$ml.pop.size
	
	s <- ecm$sample.size
	
	mpop <- ecm$max.pop
	
	get.Lik.Ratio <- function(n, s, max.pop){ # fit each data set to both ECM and TIRM
	
	dd <- simEcm(n, s)

	l.ecm <- suppressWarnings(fitEcm(dd, max.pop)$likelihood)
	
	l.tirm <- suppressWarnings(fitTirm(dd, max.pop)$likelihood)
	
	l.ratio <- l.tirm - l.ecm
	
	return(l.ratio)
	}
	
	tstat <- sapply(1:bootstraps, function(x) get.Lik.Ratio(n, s, mpop))
	
	p.value <- (length(which(tstat >= LR) + 1)) / (bootstraps + 1)
	
	return(list(LR=LR, p.value=p.value))
	
}
