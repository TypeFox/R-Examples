`snp2bin` <-
function(mat,domrec=TRUE,refAA=FALSE,snp.in.col=TRUE,monomorph=0){
	if(!is.data.frame(mat) & !is.matrix(mat))
		stop("mat must be either a matrix or a data frame.")
	mat <- as.matrix(mat)
	if(snp.in.col)
		mat<-t(mat)
	minval<-min(mat,na.rm=TRUE)
	if(!minval%in%c(0,1,"AA"))
		stop("The values of the SNPs must be either 0, 1, 2,\n or 1, 2, 3, or 'AA', 'AB', 'BB'.")
	if(minval==1 && !all(mat%in%c(1:3,NA)))
		stop("Some of the values in mat are neither 1 nor 2 nor 3.")
	if(minval==0){
		if(!all(mat%in%c(0:2,NA)))
			stop("Some of the values in mat are neither 0 nor 1 nor 2.")
		mat<-mat+1
	}
	mat.info<-if(domrec) data.frame("SNP_1"=c(0,1,1),"SNP_2"=c(0,0,1))
		else data.frame("SNP_1"=c(1,0,0),"SNP_2"=c(0,1,0),"SNP_3"=c(0,0,1))
	if(minval=="AA"){
		if(!domrec)
			refAA<-TRUE
		else
			mat.info<-data.frame(mat.info, "Assumed Genotype"=c("Homozygous Reference", 
        			"Heterozygous","Homozygous Variant"),check.names=FALSE)
		if(refAA)
			mat.info<-data.frame(SNP=c("AA","AB","BB"),mat.info)
		mat<-recodeAffySNP(mat,refAA=refAA)
	}
	else
		mat.info<-data.frame(SNP=minval+(0:2), mat.info, "Assumed Genotype"=c("Homozygous Reference", 
        			"Heterozygous","Homozygous Variant"),check.names=FALSE)
	if(is.null(rownames(mat))){
		rownames(mat)<-paste("SNP",1:nrow(mat),sep="")
		warning("Since mat has no ",ifelse(snp.in.col,"column","row")," names, generic ones are added.")
	}
	if(domrec)
		mat<-getBinDomRec(mat)
	else
		mat<-getBinGenotype(mat)
	rs<-rowSums(mat,na.rm=TRUE)
	ids<-which(rs<=monomorph | rs>=rowSums(!is.na(mat))-monomorph)
	if(length(ids)>0){
		mat<-mat[-ids,]
		warning(length(ids)," monomorphic dummy variables are removed.")
	}
	if(snp.in.col)
		mat<-t(mat)
	cat("SNPs are coded as follows:\n")
    	print(mat.info)
	mat
}

