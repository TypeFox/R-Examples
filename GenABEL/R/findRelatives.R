#' guesses relations between individuals 
#' 
#' This function guesses relationships (expressed as estimated 
#' number meiotic connection(s)) using genomic data. Compared 
#' to guessing relations from genomic kinship matrix, this 
#' procedure offers several enhancements: 
#' 
#' (1) by use of IBD/IBS 
#' 3-state space, it allows to distinguish between some pairs, 
#' which have the same kinship (e.g. parent-offspring from 
#' brother-sister; uncle-nephew from grandparent-grandchild, etc.)  
#' 
#' (2) it reports likelihood, allowing for more rigorous inferences
#' 
#' If 'gtdata' are provided as a matrix (or 'databel' matrix), 
#' genotypes should be coded as 0, 1, or 2; each SNP corresponds 
#' to a column and each ID is a row. 'q' corresponds to the frequency 
#' of 'effect' (aka 'coded') allele, which is also equivalent to the 
#' mean(SNP)/2.0 provided coding is correct.
#' 
#' 'nmeivec' is a sequence of integers, e.g. c(1,2) will test for 
#' parent-offspring pairs and pairs separated by two meioses (sibs, 
#' grandparent-grandchild, etc.). If 'nmeivec' does not contain '0' 
#' as its first element, it will be automatically added (testing for 
#' twins). Also, nmeivec will be updated with 
#' c(nmeivec,max(nmeivec)+1,100) to allow for thesting of 
#' testing vs. 'null' (unrelated, 100) and 'most distant'
#' specified by user (max(nmeivec)). 
#' 
#' While one may be interested to test only a sub-set of the data
#' for relationships, 
#' it is recommended to provide 'q' estimated using all data 
#' available (see example). 
#' 
#' 'gkinCutOff' allows use of genomic kinship matrix 
#' (computed internally) to pre-screen pairs to be tested. 
#' Use of this option with value '-1' is recommended: 
#' in this case threshod is set to 0.5^(max(nmeivec)+2).
#' If not NULL, only pairs passing gkinCutOff are tested with 
#' the likelihood procedure.   
#' 
#' After likelihood estimation, inference on relationship is made.
#' Releationship (in terms of number of meioses) is 'guessed'  
#' if odds of likelihoods under the meiotic distance providing max likelihood 
#' and under the 'null' (maximal meiotic distance tested + 100)  
#' is greater than 'OddsVsNull' parameter AND odds max-lik vs. 
#' the next-best meiotic distance   
#' is greater than 'OddsVsNextBest' parameter.
#' 
#' @param gtdata genotypic data, either 'gwaa.data' or 'snp.data' 
#' class, or matrix or 'databel' matrix (see details for format).
#' @param nmeivec vector providing the degree of relationship 
#' to be tested (1: parent-offspring; 2: sibs, grandparent-grandchild; 
#' etc.). 
#' @param q vector of effect allele frequencis for the data. 
#' @param epsilon genotyping error rate
#' @param quiet if TRUE, screen outputs supressed
#' @param OddsVsNull threshold used in relationships inferences 
#' (see details)
#' @param OddsVsNextBest  threshold used in relationships inferences 
#' (see details)
#' @param twoWayPenalty penalty on likelihoods resulting from models 
#' assuming two meiotic pathways
#' @param doTwoWay or not
#' @param vsIDs specific IDs to be tested vs others
#' @param gkinCutOff if not null, sets a threshold used 
#' to pre-screen pairs before guessing relations. If value 
#' < 0 provided, procedure sets threshold automatically (recommended)
#' @param kinshipMatrix (genomic) kinship matrix (used if gkinCutOff!=NULL) 
#' 
#' @return 
#' A list with elements
#' call -- details of the call;
#' profile -- table detailing likelihood for all pairs 
#' tested; 
#' estimatedNmeioses -- nids x nids matrix containing 
#' maximum likelihood estimate of  
#' meiotic distance for all pairs of individuals
#' guess -- same as estimatedNmeioses, but all estimates not 
#' passing inference criteria (OddsVsNull, OddsVsNextBest) are NAed;
#' compressedGuess -- same as above, but removing cows and cols 
#' with missing-only elemnts
#' 
#' @keywords htest 
#' 
#' @examples 
#' data(ge03d2.clean)
#' df <- ge03d2.clean[,autosomal(ge03d2.clean)]
#' df <- df[,sort(sample(1:nsnps(df),1000))]
#' eaf <- summary(gtdata(df))$"Q.2"
#' ### donotrun
#' \dontrun{
#' relInfo <- findRelatives(df[27:30,],q=eaf)
#' relInfo
#' # look only for 1st and 2nd degree relatives
#' relInfo1 <- findRelatives(df[27:30],q=eaf,gkinCutOff=-1,nmeivec=c(1,2,3))
#' relInfo1
#' relInfoVS <- findRelatives(df[27:30,],q=eaf,nmeivec=c(1:6),vsIDs=idnames(df[27:30,])[1:2])
#' relInfoVS
#' }
#' ### end norun
#' 
findRelatives <- function(gtdata,nmeivec=c(1:2),q=NULL,epsilon=0.01, 
		quiet=FALSE,OddsVsNull=1000,OddsVsNextBest=100,
		twoWayPenalty=log(10),doTwoWay=TRUE,vsIDs=NULL,
		gkinCutOff=NULL,kinshipMatrix=NULL) {
# sanity checks
	nmeivec <- sort(nmeivec)
	if (nmeivec[1]<0) stop("negative elements in nmeivec")
	if (nmeivec[1] != 0) nmeivec <- c(0,nmeivec)
	nmeivec <- c(nmeivec,nmeivec[length(nmeivec)]+1,nmeivec[length(nmeivec)]+100)
	if (class(gtdata) != "matrix" && class(gtdata) != "gwaa.data" && 
			class(gtdata) != "snp.data" && class(gtdata) != "databel")
		stop("gtdata should be gwaa.data, snp.data, databel or a matrix")
	if (class(gtdata) == "gwaa.data") gtdata <- gtdata(gtdata)
# get ID names
	if (!is.null(rownames(gtdata)))
		idnam <- rownames(gtdata)
	else
		idnam <- as.character(1:dim(gtdata)[1])
# get Q if not provided
	if (is.null(q)) {
		if (class(gtdata) =="snp.data") q <- summary(gtdata)$"Q.2"
		else if (class(gtdata) == "matrix") q <- colMeans(gtdata,na.rm=T)/2.
		else if (class(gtdata) == "databel") {
			q <- rep(NA,dim(gtdata)[2])
			for (i in 1:dim(gtdata)[2]) q[i] <- mean(gtdata[,1],na.rm=T)/2.
		}
		warning('"q" is not specified; inferring from data')
	}
# set vsIDs if used
	if (!is.null(vsIDs)) {
		if (is.character(vsIDs)) {
			if (!all(vsIDs %in% idnam)) stop("can not find some of vsIDs")
			vsIDs <- which(idnam %in% vsIDs)
			if (length(vsIDs)>1) {ma <- max(vsIDs); mi <- min(vsIDs)} 
			else ma <- mi <- vsIDs
			if (ma>dim(gtdata)[1] || mi<1) stop("vsIDs out of range")
		}
		notVsIDs <- c(1:dim(gtdata)[1])[which(!(c(1:dim(gtdata)[1]) %in% vsIDs))]
	}
# get gkin if used
	lengthOut <- dim(gtdata)[1]
	if (is.null(vsIDs)) lengthOut <- lengthOut*(lengthOut-1)/2
	else lengthOut <- (lengthOut-length(vsIDs))*length(vsIDs)
	if (!is.null(gkinCutOff)) {
		if (is.null(kinshipMatrix)) {
			if (is(gtdata,"snp.data")) {
				if (is.null(q)) {
					kinshipMatrix <- ibs(gtdata,weight="freq")
				} else {
					kinshipMatrix <- ibs(gtdata,weight="freq",snpfreq=q)
				}
			}
			else stop("gkinCutOff can only be used with gtdata of gwaa.data or snp.data class")
		} else {
			kinshipMatrix <- kinshipMatrix[idnam,idnam]
		}
		if (gkinCutOff<(-0.5)) {
			gkinCutOff <- .5^(max(nmeivec[1:(length(nmeivec)-1)])+1)
			if (!quiet) cat("*** gkinCutOff set to ",gkinCutOff,"\n")
		}
		if (is.null(vsIDs)) {
			testUs <- as.vector(kinshipMatrix[lower.tri(kinshipMatrix)])
			testUs <- (testUs>gkinCutOff)
		} else {
#			testUs <- c()
#			for (jj in vsIDs) {
#				if (jj>1) 
#				testUs <- c(testUs,as.vector((kinshipMatrix[jj,])[1:(jj-1)]),
#						as.vector((kinshipMatrix[,jj])[c((jj+1):dim(kinshipMatrix)[1])]))
#			else 
#				testUs <- c(testUs,
#						as.vector((kinshipMatrix[,jj])[c((jj+1):dim(kinshipMatrix)[1])]))
#			testUs <- 
#				print(jj)
#				print(testUs)
#			}
			kinshipMatrix[upper.tri(kinshipMatrix)] <- t(kinshipMatrix)[upper.tri(kinshipMatrix)]
			testUs <- kinshipMatrix[vsIDs,notVsIDs]
#			print(testUs)
			testUs <- (testUs>gkinCutOff)
		}
		if (length(testUs) != lengthOut) 
			stop(paste("length(testUs) != lengthOut:",
							length(testUs),";",lengthOut))
	} else {
		testUs <- rep(TRUE,lengthOut)
	}
	if (all(testUs==FALSE)) stop("no test to be performed according to specified criteria")
# prepare 2-way meiotic table to iterate over
	meiTab <- list()
	#print(nmeivec)
	maxNMeiNotNullNotBoundary <- nmeivec[(length(nmeivec)-2)]
	for (i in nmeivec[1:(length(nmeivec)-1)]) {
		fromTo <- c(max(2,i):maxNMeiNotNullNotBoundary)
		if (doTwoWay && i>0 && i <= maxNMeiNotNullNotBoundary) {
			series <- matrix(c(rep(i,length(fromTo)),fromTo),ncol=2)
			if (dim(series)[1]>=1)
				for (j in 1:dim(series)[1]) {
					#print(c(i,j))
					cname <- paste(series[j,1],series[j,2],sep="+")
					meiTab[[cname]] <- series[j,]
				}
		}
		meiTab[[as.character(i)]] <- i
	}
	meiTab[[ as.character(nmeivec[length(nmeivec)]) ]] <- nmeivec[length(nmeivec)]
	#print(meiTab)
	#print(names(meiTab))
	#stop()
# iterate over IDs
	outI <- 1
	out <- matrix(ncol=(length(meiTab)),nrow=lengthOut) 
	outN1 <- matrix(ncol=1,nrow=lengthOut)
	outN2 <- matrix(ncol=1,nrow=lengthOut)
	if (!quiet) pb <- txtProgressBar(style = 3)
	if (is.null(vsIDs)) {id1list <- c(1:(dim(gtdata)[1]-1))}
	else { id1list <- vsIDs }
	testsDone <- 0
	for (id1 in id1list) {
		name1 <- idnam[id1]
		if (class(gtdata) == "snp.data")
			gt1 <- as.numeric(gtdata[id1,])
		else 
			gt1 <- gtdata[id1,]
		bgt1 <- blurGenotype(gt1,q=q,epsilon=epsilon)
		if (is.null(vsIDs)) id2list <- (id1+1):dim(gtdata)[1]
		else id2list <- notVsIDs #c(c(1:(id1-1)),c((id1+1):dim(gtdata)[1]))
		for (id2 in id2list) {
			name2 <- idnam[id2]
			#print(id1,id2)
			#print(c(outI,testUs))
			if (testUs[outI]) { 
				if (class(gtdata) == "snp.data")
					gt2 <- as.numeric(gtdata[id2,])
				else 
					gt2 <- gtdata[id2,]
				bgt2 <- blurGenotype(gt2,q=q,epsilon=epsilon)
				acuL <- rep(NA,length(nmeivec)); 
				k <- 1
				for (mei in meiTab) {
					# AAA
					TM <- makeTransitionMatrix(q,nmeioses=mei)
					acuL[k] <- getLogLikelihoodGivenRelation(bgt1,bgt2,TM)$logLik
					if (length(mei)>1) acuL[k] <- acuL[k] - twoWayPenalty
					k <- k+1
				}
				testsDone <- testsDone + 1
				if (!quiet) setTxtProgressBar(pb, testsDone/sum(testUs))
			} else {
				acuL <- c(1:length(meiTab))
			}
			outN1[outI,] <- name1
			outN2[outI,] <- name2
			#cat("id1list",id1list,"\n")
			#cat("id2list",id2list,"\n")
			#cat("\n",acuL,"\n")
			out[outI,] <- c(acuL[1:(length(acuL))]-acuL[length(acuL)])
			#print(out[outI,])
			outI <- outI+1
		}
	}
	if (!quiet) cat("\n")
	#print(out)
	meiMtmp <- apply(out,MARGIN=1,FUN=function(x){return(which.max(x))})
	#print(meiMtmp)
	if (is.null(vsIDs)) {
		meiM <- matrix(ncol=dim(gtdata)[1],nrow=dim(gtdata)[1])
		diag(meiM) <- 0
		meiM[lower.tri(meiM)] <- names(meiTab)[meiMtmp]
		meiM[upper.tri(meiM)] <- t(meiM)[upper.tri(meiM)]
		colnames(meiM) <- rownames(meiM) <- idnam
	} else { 
		meiM <- matrix(ncol=length(notVsIDs),nrow=length(vsIDs))
		meiM[] <- names(meiTab)[meiMtmp]
		rownames(meiM) <- idnam[vsIDs]
		colnames(meiM) <- idnam[notVsIDs]
	}
	#print(meiM)
	firstBestLik <- apply(out,MARGIN=1,FUN=function(x){return(max(x))})
	which_firstBestLik <- apply(out,MARGIN=1,FUN=function(x){return(which.max(x))})
	secondBestLik <- apply(out,MARGIN=1,FUN=function(x){x[order(x,decreasing=TRUE)[2]]})
	cndX <- !(exp(firstBestLik-secondBestLik)>OddsVsNextBest 
				& exp(firstBestLik)>OddsVsNull & (which_firstBestLik != (length(meiTab)-1)))
	if (is.null(vsIDs)) {
		cnd <- diag(FALSE,dim(gtdata)[1],nrow=dim(gtdata)[1])
		cnd[lower.tri(cnd)] <- cndX
		cnd[upper.tri(cnd)] <- t(cnd)[upper.tri(cnd)]
	} else {
		cnd <- matrix(ncol=length(notVsIDs),nrow=length(vsIDs))
		cnd[] <- cndX
	}
	guess <- meiM
	guess[which(cnd==1)] <- NA  
	out <- data.frame(outN1,outN2,out,stringsAsFactors=FALSE)
	names(out) <- c("id1","id2",names(meiTab))
	tmp <- guess; diag(tmp) <- NA; 
	todrop <- apply(tmp,MARGIN=1,FUN=function(x){return(all(is.na(x)))})
	compressedGuess <- guess[!todrop,!todrop] 
	finalOut <- list(call=match.call(),profile=out,estimatedNmeioses=meiM,
			firstBestLik=firstBestLik,
			secondBestLik=secondBestLik,guess=guess,compressedGuess=compressedGuess)
	finalOut
}
