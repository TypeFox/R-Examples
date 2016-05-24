function.ROC01 <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
		
	distance <- (measures.acc$Sp[,1]-1)^2+(measures.acc$Se[,1]-1)^2
	cROC01 <- measures.acc$cutoffs[which(round(distance,10) == round(min(distance,na.rm=TRUE),10))]
	optimal.distance <- min(distance,na.rm=TRUE)	

	optimal.cutoff <- obtain.optimal.measures(cROC01, measures.acc)

	res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = distance, optimal.criterion = optimal.distance)
	res
}
