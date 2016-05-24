summary.credibleshiftset <- function(object, ...) {

	conf <- round(100*object$set.limit);
	cat('\n', conf, '% credible set of rate shift configurations sampled with BAMM');
	cat('\n\nDistinct shift configurations in credible set: ', object$number.distinct);
	
	mm <- min(c(9, object$number.distinct));
	if (object$number.distinct > 9){
		omitted <- object$number.distinct - 9;
	}
	
	xvec <- numeric(object$number.distinct);
	dd <- data.frame(rank=1:object$number.distinct, probability = xvec, cumulative=xvec, Core_shifts=xvec);
	for (i in 1:nrow(dd)){
		dd$probability[i] <- object$frequency[i];
		dd$cumulative[i] <- object$cumulative[i];
		dd$Core_shifts[i] <- length(object$shiftnodes[[i]]);
	}
	
	
	cat('\n\nFrequency of', mm, 'shift configurations with highest posterior probability:\n\n\n');
	write.table(format(t(names(dd)),justify="centre",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);
	write.table(format(dd[1:mm,], justify="none", width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);
	
#	cat('\n\tRank\t\tProb\tCumulative\t\tCoreShifts\n')
#	for (i in 1:mm){
#		cat('\t\t', format(i, justify='left'), '\t\t', format(round(object$frequency[i], 3), nsmall=3), '\t');
#		cat(round(object$cumulative[i], 3), '\t\t\t', length(object$shiftnodes[[i]]), '\n');
#	}	
	
	if (object$number.distinct > 9){
		cat('\n...omitted', omitted, 'additional distinct shift configurations\n');
		cat('from the credible set. You can access the full set from your \n');
		cat('credibleshiftset object\n');
	}
	
	cat('\n')
	

	
	invisible(dd);
}


