### --- Test setup ---

if(FALSE) {
	## Not really needed, but can be handy when writing tests
	library("RUnit")
	library("GenABEL")
	library("DatABEL")
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source(paste("../inst/unitTests/shared_functions.R"))
source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.impute2databel <- function()
{
	
	#library(GenABEL)
	#library("RUnit")
	unlink("*.fv?")
	
	makedose <- function(prob) {
		dose <- 2*prob[c(F,F,T)]+prob[c(F,T,F)]
		bp <- prob[c(T,F,F)]
		miss <- which((bp+dose)<0.9)
		if (length(miss) > 0 ) dose[miss] <- NA
		return(dose)
	}
	makeprob <- function(prob) {
		dose2 <- prob[c(F,F,T)]
		dose1 <- prob[c(F,T,F)]
		bp <- prob[c(T,F,F)]
		miss <- which((bp+dose1+dose2)<0.9)
		out <- matrix(c(dose2,dose1),ncol=2)
		if (length(miss) > 0 ) out[miss,] <- NA
		out
	}
	
	
	tmp0 <- read.table("TEST10x15.geno",head=F,strings=F)
	snps <- tmp0[,2]
	tmp0 <- tmp0[,c(6:dim(tmp0)[2])]
	dose <- apply(tmp0,FUN="makedose",MAR=1)
	rownames(dose)
	prb <- makeprob(tmp0[1,])
	for (i in 2:dim(tmp0)[1])
		prb <- cbind(prb,makeprob(tmp0[i,]))
	dm <- dim(prb)
	prb <- as.numeric(prb)
	prb <- matrix(prb,ncol=dm[2])
	
	tmp1 <- impute2databel(geno="TEST10x15.geno",
			sample="impute.sample5",
			out="TEST10x15_T.geno",
			makeprob=FALSE,
			old=TRUE)
	tmp1
	tmp1_m <- as(tmp1,"matrix")
	tmp1_m
	
	
	print(tmp0[1:5,1:6])
	tmp2 <- impute2databel(geno="TEST10x15.geno",
			sample="impute.sample5",
			out="TEST10x15_F.geno",
			makeprob=TRUE,
			old=FALSE)
	print(tmp2)
	tmp2_m <- as(tmp2,"matrix")
	
	
	tmp3 <- databel("TEST10x15_F.geno.prob")
	tmp3_m <- as(tmp3,"matrix")
	
	checkIdentical(TRUE,all(abs(dose-tmp1_m)<1e-7,na.rm=TRUE))
	checkIdentical(TRUE,all(abs(dose-tmp2_m)<1e-7,na.rm=TRUE))
	checkIdentical(TRUE,all(abs(prb-tmp3_m)<1e-7,na.rm=TRUE))
	checkIdentical(TRUE,all(abs(tmp1_m-tmp2_m)<1e-8,na.rm=TRUE))
	checkEqualsNumeric(dose,tmp1_m,checkNames=FALSE,tolerance = 5*(.Machine$double.eps^0.5))
	checkEqualsNumeric(dose,tmp2_m,checkNames=FALSE,tolerance = 5*(.Machine$double.eps^0.5))
	checkEqualsNumeric(prb,tmp3_m,checkNames=FALSE,tolerance = 5*(.Machine$double.eps^0.5))
	checkIdentical(tmp1_m,tmp2_m)
	checkEqualsNumeric(tmp1_m,tmp2_m)
	
	smpl <- read.table("impute.sample5",head=F,skip=2,strings=F)[,1]
	rownames(dose) <- smpl
	rownames(prb) <- smpl
	
	checkEquals(dose,tmp1_m,tolerance=4*sqrt(.Machine$double.eps))
	checkEquals(dose,tmp2_m,tolerance=4*sqrt(.Machine$double.eps))
	checkEquals(prb,tmp3_m,tolerance=4*sqrt(.Machine$double.eps))
	checkIdentical(tmp1_m,tmp2_m)
	checkEqualsNumeric(tmp1_m,tmp2_m)
	checkIdentical(rownames(tmp1),smpl)
	checkIdentical(rownames(tmp2),smpl)
	checkIdentical(rownames(tmp3),smpl)
	checkIdentical(dimnames(tmp1),dimnames(tmp1_m))
	checkIdentical(dimnames(tmp2),dimnames(tmp2_m))
	checkIdentical(dimnames(tmp3),dimnames(tmp3_m))
	checkIdentical(get_dimnames(tmp1),list(smpl,snps))
	checkIdentical(get_dimnames(tmp2),list(smpl,snps))
	snps2 <- rep(0,2*length(snps))
	snps2[c(T,F)] <- paste(snps,"_11",sep="")
	snps2[c(F,T)] <- paste(snps,"_01",sep="")
	checkIdentical(get_dimnames(tmp3),list(smpl,snps2))
	
	rm(list=ls());gc()
	
	unlink("*.fv?")
	
}

# fixing CRAN NOTE
#
test.SNPquoted <- function() {	
	#library("RUnit")
	#library("GenABEL")
	require(DatABEL)
# generate some random data in "prob" format
	Nids <- 10
	Nsnps <- 3
	q <- runif(Nsnps,min=0.2,max=0.8)
	PROB0 <- apply(matrix(q,ncol=1),MARGIN=1,FUN=function(x){
				cBB <- rnorm(Nids,mean=x^2,sd=0.01)
				cAB <- rnorm(Nids,mean=2*x*(1-x),sd=0.01)
				cAA <- 1-cBB-cAB
				cAA <- matrix(cAA,ncol=1)
				return(c(cAA,cAB,cBB))
			})
	PROB <- matrix(NA,ncol=3*Nsnps,nrow=Nids)
	for (i in 1:Nsnps) {
		for (j in 1:3) {
			PROB[,(i-1)*3+j] <- PROB0[c(((j-1)*Nids+1):(j*Nids)),i]
		}
	}
	daProb <- matrix2databel(t(PROB),filename="temp_daProb")
	
# conversion functions
	makedose <- function(prob) {
		dose <- 2*prob[c(F,F,T)]+prob[c(F,T,F)]
		bp <- prob[c(T,F,F)]
		miss <- which(abs(bp)<1e-16 & abs(dose)<1e-16)
		if (length(miss) > 0 ) dose[miss] <- NA
		return(dose)
	}
	pfun <- function(a) return(a)

	dosefile <- apply2dfo(dfodata=daProb, anFUN = "makedose", 
			MAR = 2, procFUN = "pfun",prob=get("SNP"),
			outclass="databel",
			outfile="temp_myTestDose",
			type="DOUBLE",transpose=FALSE)
	
# check equality
	DOSE <- matrix(NA,ncol=Nids,nrow=Nsnps)
	for (i in 1:Nids) for (j in 1:Nsnps) {
			DOSE[j,i] <- 2*PROB[i,j*3]+PROB[i,j*3-1]
		}
	
	checkEqualsNumeric(DOSE,as(dosefile,"matrix"))
	
	rm(dosefile); gc()
	unlink("temp_myTestDose.fv?")
	
	rm(daProb); gc()
	unlink("temp_daProb.fv?")
}
