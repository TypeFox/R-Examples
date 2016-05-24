"HWE.show" <-
function(data,idsubset=c(1:data@gtdata@nids),snpsubset =c(1:data@gtdata@nsnps)) {
	if (is(data,"gwaa.data")) {
        	if (missing(idsubset))  idsubset  <- c(1:data@gtdata@nids)
        	if (missing(snpsubset)) snpsubset <- c(1:data@gtdata@nsnps)
		n1data <- data@gtdata[idsubset,snpsubset]
	} else if (is(data,"snp.data")) {
        	if (missing(idsubset))  idsubset  <- c(1:data@nids)
        	if (missing(snpsubset)) snpsubset <- c(1:data@nsnps)
		n1data <- data
	} else {
		stop("data should be of type \"gwaa.data-class\" or \"snp.data-class\"")
	}
	for (i in 1:length(snpsubset)) {
		if (n1data[,snpsubset[i]]@chromosome == "X") {
			nidsub <- (n1data[,snpsubset[i]]@male==0)
			cat ("This is X-chromsome marker. Female data shown...\n")
		} else
			nidsub <- n1data[,snpsubset[i]]@idnames 
		gt <- as.double(n1data[nidsub,snpsubset[i]])
		gtc<- as.character(n1data[nidsub,snpsubset[i]])
		out <- matrix(NA,nrow=3,ncol=3)
		out[1,1] <- sum(1*(gt==0),na.rm=TRUE)
		out[1,2] <- sum(1*(gt==1),na.rm=TRUE)
		out[1,3] <- sum(1*(gt==2),na.rm=TRUE)
		N <- sum(out[1,])
		sumdat <- summary(n1data[,i])
		f2 <- sumdat[1,"Q.2"]
		f1 <- 1 - f2
		colnames(out) <- c("A/A","A/B","B/B")
		out[2,] <- c(f1*f1*N,2*f1*f2*N,f2*f2*N)
		if (!any(out[2,]<=0))
			out[3,] <- (out[1,]-out[2,])^2/out[2,]
		else
			out[3,] <- (out[1,]-out[2,])^2
		rownames(out) <- c("observed","expected","chi2+")
		chi2 <- sum(out[3,])	
		pv <- pchisq(chi2,1,lower.tail=FALSE)
		sco <- cat("HWE summary for",snpsubset[i],":\n")
		sco <- print(out)
		sco <- cat("Chi2 =",chi2,"; P =",pv)
		pvex <- sumdat[1,"Pexact"]
		sco <- cat("; exact P =",pvex,"\n\n")
		sco
	}
}

