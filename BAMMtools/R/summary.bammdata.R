summary.bammdata = function(object, display=10, print=T, ...) {

	fev <- sapply(object$eventData, nrow);
	xx <- table(fev - 1);
	shifts <- as.numeric(names(xx));
	xx <- as.numeric(xx);
	df <- data.frame(shifts = shifts, prob = signif(xx / sum(xx), digits= 3));

	if (print){
		
		cat("\nAnalyzed", length(object$eventData), "posterior samples\n");
	 
		cat("Shift posterior distribution:\n\n");
		disp <- data.frame(cbind(df$shifts),signif(df$prob,2));
		if (nrow(disp) <= display) {	
			write.table(format(disp,justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);	
		} else {
			if(nrow(disp)-display == 1) {
				write.table(format(disp[1:display,],justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);
				cat("... omitted 1 row\n\n");
			}else {
				wr <- which(disp[,2] == max(disp[,2]))[1];
				index <- c(max(1,floor(wr-display/2)), wr, min(ceiling(wr+display/2),nrow(disp)));
				if (index[1] > 1) {
					cat("... omitted",index[1]-1,"rows\n");
				}
				write.table(format(disp[index[1]:index[3],],justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);
				cat("... omitted", nrow(disp)-index[3]-1,"rows\n\n");
			}
		}
		cat("\nCompute credible set of shift configurations for more information:\n");
		cat("\tSee ?credibleShiftSet and ?getBestShiftConfiguration\n");		
		
	}
	
	invisible(df);
}
