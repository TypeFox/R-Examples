plot.ei_compare <-
function(x, ...) {
	
	
	diff_table <- x@data[, grep("EI_Diff", colnames(x@data))]
	diff_table <- na.omit(diff_table)
	
	cand_names <- as.character ( x@data[,1] )
	cand_names <- cand_names[seq(1,length(cand_names),2)]

	group_names <- x@groups
	
	nplots <- ncol(diff_table) # Number of panels to make
	
	
	if (nplots==1) {
		diff_table <- diff_table[,1]
	
		if (max(abs(diff_table)) < 10) {
			xlim <- c(-10,10)
		} else {
			xlim <- c ( - (ceiling ( max(abs(diff_table))) + 2), ceiling ( max(abs(diff_table))) + 2 )
		}
		plot(diff_table, 1:length(diff_table), xlim=xlim, pch=19, xlab = "RxC - EI", ylab = "Candidates",
			main = paste(group_names, "EI Diff Comparison", sep=":"), yaxt="n", ...)
		axis(2, at = 1:length(diff_table), labels = cand_names, las=1, ...) # need to adjust sideways

		abline(v=0, lty=2, lwd=3, col="grey")
		abline(v=1, lty=2, lwd=1, col="grey")
		abline(v=-1, lty=2, lwd=1, col="grey")
		abline(v=5, lty=2, lwd=1, col="grey")
		abline(v=-5, lty=2, lwd=1, col="grey")
	
	} else if ( nplots == 2) {
	
		par(mfrow=c(2,1))
	
		for (i in 1:nplots) {
	
			if (max(abs(diff_table[,i])) < 10) {
				xlim <- c(-10,10)
			} else {
				xlim <- c ( - (ceiling ( max(abs(diff_table[,i]))) + 2), ceiling ( max(abs(diff_table[,i]))) + 2 )
			}
			
			plot(diff_table[,i], 1:length(diff_table[,i]), xlim=xlim, pch=19, xlab = "RxC - EI", ylab = "Candidates",
				main = paste(group_names[i], " EI Diff Comparison",sep=":"), yaxt="n", ...)
			axis(2, at = 1:length(diff_table[,i]), labels = cand_names, las=1, ...) # need to adjust sideways

			abline(v=0, lty=2, lwd=3, col="grey")
			abline(v=1, lty=2, lwd=1, col="grey")
			abline(v=-1, lty=2, lwd=1, col="grey")
			abline(v=5, lty=2, lwd=1, col="grey")
			abline(v=-5, lty=2, lwd=1, col="grey")
		}

	} else if ( nplots == 3 | nplots == 4) {
	
		par(mfrow=c(2,2))
	
		for (i in 1:nplots) {
	
			if (max(abs(diff_table[,i])) < 10) {
				xlim <- c(-10,10)
			} else {
				xlim <- c ( - (ceiling ( max(abs(diff_table[,i]))) + 2), ceiling ( max(abs(diff_table[,i]))) + 2 )
			}
			
			plot(diff_table[,i], 1:length(diff_table[,i]), xlim=xlim, pch=19, xlab = "RxC - EI", ylab = "Candidates",
				main =  paste(group_names[i], " EI Diff Comparison",sep=":"), yaxt="n", ...)
			axis(2, at = 1:length(diff_table[,i]), labels = cand_names, las=1, ...) # need to adjust sideways

			abline(v=0, lty=2, lwd=3, col="grey")
			abline(v=1, lty=2, lwd=1, col="grey")
			abline(v=-1, lty=2, lwd=1, col="grey")
			abline(v=5, lty=2, lwd=1, col="grey")
			abline(v=-5, lty=2, lwd=1, col="grey")
		} # Close nplots for loop

	} # Close else if loop
	
}
