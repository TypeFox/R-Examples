print.linkdat <- function(x, ..., markers) {
	if (missing(markers)) marker.nos = seq_len(min(x$nMark, 5)) 
	else {
		if(length(markers) > 0 && max(markers) > x$nMark) stop("Nonexisting marker(s) indicated") 
		marker.nos = markers
	}
	datafr = as.data.frame(x, markers=marker.nos, sep="/", missing="-", singleCol=TRUE)
   print(datafr, ...)
	if (missing(markers) && x$nMark > 5)
		cat("\nOnly first 5 markers are shown. Use option 'markers=' to print specified markers.\n")
   invisible(datafr)
}

print.linkdat.model <- function(x, ...) {
	model = x
	switch(model$chrom, 
		AUTOSOMAL = cat("Autosomal inheritance with penetrances: (f0, f1, f2) =", paste('(',paste(model$penetrances, collapse=", "),')',sep=""), "\n"),
		X = cat("X-linked inheritance with penetrances:\n\tMales: (f0, f1) =", paste('(', paste(model$penetrances$male, collapse=", "),')',sep=""), 
				"\n\tFemales: (f0, f1, f2) =", paste('(', paste(model$penetrances$female, collapse=", "),')',sep=""), "\n")
	)
	cat("Disease allele frequency:", model$dfreq,"\n")
}

summary.linkdat <- function(object, ...) {
	x <- object
	cat("Pedigree:\n---------\n")
	cat(x$nInd,"individuals\n")
	cat(length(x$founders),"founders,", length(x$nonfounders),"nonfounders; bit size =", 2*length(x$nonfounders)-length(x$founders), "\n")
	if((ant<-length(x$subnucs)) > 0) cat(ant,"nuclear", ifelse(ant==1, "subfamily","subfamilies"),"\n")
	aff = x$pedigree[, 'AFF']
	if(all(aff==1)) cat("No pedigree members affected by disease\n")
	else cat(sum(aff==2), "affected by disease,", sum(aff==1), "unaffected,", sum(aff==0), "with unknown affection status\n") 
	
	cat("\nMarker data:\n------------\n", x$nMark, ifelse(x$nMark==1, " marker ", " markers "), "in total\n", sep="")
	if (x$nMark > 0) {
		miss = which(rowSums(m <- do.call(cbind, x$markerdata)) == 0)
		cat(length(miss), "individuals with no available genotypes")
		if(length(miss)>0) {
			cat(":", paste(x$orig.ids[miss], collapse=", "), "\n")
			cat(round(sum(m[-miss,]==0)/length(m[-miss,])*100, 2),"% missing alleles (excluding ungenotyped individuals)\n")
		} else 
			cat("\n", sum(m==0)/length(m)*100, "% missing alleles\n", sep="")
		cat("\nChromosome distribution of markers:\n")
		chrtbl = table(sapply(x$markerdata, attr, 'chrom'), useNA="ifany")
		names(chrtbl)[is.na(names(chrtbl))] = "unknown"
		for(i in seq_along(chrtbl)) 
			cat(" chromosome ", names(chrtbl)[i], ": ", chrtbl[i], ifelse(chrtbl[i]==1, " marker\n", " markers\n"), sep="")
		
		cat("\nAllele number distribution:\n")
		ntbl = table(unlist(lapply(x$markerdata, attr, 'nalleles')))
		for(i in seq_along(ntbl)) 
			cat(" ", names(ntbl)[i], " alleles", ": ", ntbl[i], ifelse(ntbl[i]==1, " marker\n", " markers\n"), sep="")
	}

	
	cat("\nModel parameters:\n-----------------\n")
	if(is.null(x$model)) cat("No model parameters set\n")
	else print(x$model)
}
