summary.analytics_virgin <- function(object, ...) {	
	l <- object@label_summary
	
	total <- colSums(l)[1]
	prp_cc <- round(l[,1]/total*100, digits=2)
	prp_pc <- round(l[,2]/total*100, digits=2)
	
	out <- cbind(l,PROP_CONSENSUS_CODED=prp_cc,PROP_PROBABILITY_CODED=prp_pc)

	cat("SUMMARY OF LABEL CLASSIFICATION\n\n")
	print(out)
}