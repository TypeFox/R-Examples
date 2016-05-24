"descriptives.scan" <-
function(data,file,top=10,sortby="P1df",sep="\t") {
	if (!is(data,"scan.gwaa")) stop("data argument must be of scan.gwaa-class")

	sbargs <- c("no","P1df","P2df","Pc1df")

	if (!(sortby %in% sbargs)) {cat("sortby argument must be one of the following:",sbargs);stop("");}

	#data$P1df[is.na(data[,"P1df"])] <- 9.99
	#data$P2df[is.na(data[,"P2df"])] <- 9.99
	#data$Pc1df[is.na(data[,"Pc1df"])] <- 9.99

	if (sortby=="P1df") ix <- order(data[,"P1df"])[1:top] #,index.return=TRUE)$ix[1:top]
	else if (sortby=="P2df") ix <- order(data[,"P2df"])[1:top] #,index.return=TRUE)$ix[1:top]
	else if (sortby=="Pc1df") ix <- order(data[,"Pc1df"])[1:top] #,index.return=TRUE)$ix[1:top]
	else if (sortby=="no") {ix <- c(1:min(top,dim(data)[1]))}
	else {cat("sortby argument must be one of the following:",sbargs);stop("");}

	out <- results(data)[ix,]
	#out[c(1:top),1] <- data[ix,"N"]
	#out[c(1:top),2] <- data[ix,"effB"]
	#out[c(1:top),3] <- data[ix,"P1df"]
	#out[c(1:top),4] <- data[ix,"Pc1df"]
	#out[c(1:top),5] <- data[ix,"effAB"]
	#out[c(1:top),6] <- data[ix,"effBB"]
	#out[c(1:top),7] <- data[ix,"P2df"]
	#wonder if that will work any better?!
	##out <- round(out,digits=digits)
	#out <- data.frame(out)
	#out <- cbind(map(data)[ix],out)
	#out <- cbind(chromosome(data)[ix],out)
	#rownames(out) <- snpnames(data)[ix]
	#colnames(out) <- c("Chromosome","Position","N","effB","P1df","Pc1df","effAB","effBB","P2df")
	if (!missing(file)) {
		cat(sep,file=file,sep="")
		cat(colnames(out),file=file,sep=sep,append=TRUE)
		cat("\n",file=file,sep="",append=TRUE)
		write.table(out,file=file,sep=sep,append=T,col.names=FALSE)
	}
	cat("Summary for top",top,"results, sorted by",sortby,"\n")
	out
}
