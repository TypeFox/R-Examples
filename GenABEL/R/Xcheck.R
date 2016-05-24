"Xcheck" <-
function(data,Pgte=0.01,Pssw=0.01,Pmsw=0.01,odds=1000,tabonly=FALSE,Fmale=0.8,Ffemale=0.2) {
	if (!is(data,"snp.data")) stop("data argument should be of snp.data-class")
	if (any(data@chromosome != "X")) stop("All markers should be X-linked")
	male <- (data@male==1)
	out <- list()
	out$xerr <- 0
	q <- summary(data)[,"Q.2"]
	if (sum(male)) {
		xdat <- as.numeric(data[male,])
		if (any(xdat==1,na.rm=T)) {
			out$xerr <- 1
			out$Xerrtab <- crnames(dimnames(xdat),which(xdat==1))
			colnames(out$Xerrtab) <- c("ID","SNP")
		}
	}
	if (tabonly) return(out)
	xdat <- as.numeric(data)
	xdat <- t(xdat)
	ll.female <- log(q*q*(1-Pgte)*(xdat==2)+(2*q*(1-q)*(1-Pgte)+Pgte)*(xdat==1)+(1-q)*(1-q)*(1-Pgte)*(xdat==0))
	ll.male <- log(q*(1-Pgte)*(xdat==2)+Pgte*(xdat==1)+(1-q)*(1-Pgte)*(xdat==0))
	ll.sex <- ll.female
	ll.sex[,male] <- ll.male[,male]
# find SNPs with are likely to be (pseudo)autosomal
	snpprob.0 <- apply(ll.sex,MARGIN=1,FUN=sum,na.rm=T)+log(1-Pmsw)
	snpprob.1 <- apply(ll.female,MARGIN=1,FUN=sum,na.rm=T)+log(Pmsw)
	snpODDs <- snpprob.1-snpprob.0
	out$Xmrkfail <- names(snpprob.1[snpODDs>log(odds)])
# find male which are likely to be female
	idprob.0 <- apply(ll.sex,MARGIN=2,FUN=sum,na.rm=T)+log(1-Pssw)
	idprob.1 <- apply(ll.female,MARGIN=2,FUN=sum,na.rm=T)+log(Pssw)
	idODDs <- idprob.1-idprob.0
	out$isfemale <- names(idprob.1[idODDs>log(odds)])
# find female which are likely to be male
	idprob.1 <- apply(ll.male,MARGIN=2,FUN=sum,na.rm=T)+log(Pssw)
	idODDs <- idprob.1-idprob.0
	out$ismale <- names(idprob.1[idODDs>log(odds)])
# find people with strange F
	pis <- perid.summary(data)
	out$otherSexErr <- rownames(pis)[pis$F > Ffemale & pis$F < Fmale]
# return object
	out$Xidfail <- unique(c(out$ismale,out$isfemale,out$othersexErr))
	out
}

