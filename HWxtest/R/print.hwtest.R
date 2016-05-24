# For printing the output of hwx.test
# (c) William R. Engels, 2014


#' S3 Method for printing \code{hwtest} objects

#' 
#' Prints test results (\code{hwtest}) objects depending on how much detail is provided. If histogram data are present, \code{ggplot2} is used to draw the plot by calling \code{\link{makeHistogram}}
#' 
#' @param x the results from a call to \code{\link{hwx.test}}
#' @param detail 0 for no print; 1 for P value only; 2 for all four P values; 3 to add data; 4 to add expected values
#' @param statName which statistic to use
#' @param \dots other parameters passed to \code{print}.
#' @param plotHisto Indicate whether or not to plot the histogram. Only used if \code{\link{hwx.test}} was called with \code{histobins} set to a positive value.
#' 

#' @export
print.hwtest <- 
function(x, detail=NA, statName=NA, plotHisto=TRUE, ...) {
	h <- x
	if(is.na(h$method)) return()
	if(!is.na(detail)) h$detail=detail
	if(!is.na(statName))h$statName=statName
	if(h$method=="exact" && h$tableCount < 0) {
		cat("\nCalculation timed out. You can change the time limit by increasing parameter 'safeSecs'\n")
		return();
		}
	p <- h$Pvalues;
	ob <- h$observed;
	statNames <- names(h$SE)
	comments <- character(4);
	if(ob[3]>=0) {comments[[3]] <- " (test for homozygote excess)"}
	if(ob[3]<0) {comments[[3]] <- " (test for heterozygote excess)"};
	names(comments) <- statNames;
	detail <- h$detail
	if(detail==1) {
		cat("P value (", h$statName,") = ", formatC(p[h$statName]), sep="");
		if(h$method=="monte") cat( " \u00b1 ",formatC(h$SE[h$statName], digits=5), sep="");
		cat(comments[h$statName],"\n");
	}
	if(detail >= 2) {
		cat("\n*****    Sample of ", sum(h$alleles)/2," diploids with ", length(h$alleles), " alleles", sep="")
		if(!is.null(x$sampleName)) cat(" :  ", x$sampleName)
		if(h$method=="monte") cat("\nMonte Carlo test for HW with ", h$ntrials," trials.\n", sep="");
		if(h$method=="exact") cat("\nFull enumeration of ",h$tableCount, " tables to test for HW\n", sep="");
		for(i in 1:4){
			cat("\nP value (",statNames[i],")",strtrim("     ", 6-nchar(statNames[i])),"= ", formatC(p[i], digits=6, format="f"), sep="");
			if(h$method=="monte") cat( " \u00b1 ",formatC(h$SE[i], digits=5), sep="");
			cat(comments[i]);
		}
	}
	if(detail>=3) {
		cat("\n\nObserved Test Statistics:\n");
		for(i in 1:4){cat("\n",strtrim("     ", 6-nchar(statNames[i])), statNames[i],"  :  ", formatC(ob[i]))}
	}
	if(detail>=4){
		cat("\n\nObserved Allele Counts: ", h$alleles);
		cat("\n\nObserved Genotype Counts\n");
		print(clearUpper(h$genotypes), na.print="");
	}
	if(detail>=5){
		cat("\n\nExpected Genotype Counts\n");
		ecounts <- observedX2(h$genotypes, returnExpected=T);
		rownames(ecounts) <- rownames(h$genotypes);
		colnames(ecounts) <- rownames(h$genotypes);
		print(clearUpper(ecounts), digits=3, na.print="");
	}
	cat("\n", sep="");
	if(x$histobins && plotHisto) print(makeHistogram(x))
}

