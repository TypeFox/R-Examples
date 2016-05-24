# Hinda Haned
# January 2010, Lyon

simPCR2<-function(ncells,probEx,probAlq,probPCR,cyc=28,Tdrop=2*10^7,probSperm=0.5,dip=TRUE,KH=55)
{
#standard verifications
	
	if(is.null(ncells) || !is.numeric(ncells) || ncells <0)
	{
		stop("'ncells' must be a numeric giving the number of cells")
	}
	if(cyc<1){stop("Th number of PCR cycles must at least equal 1")}
	if(Tdrop <10^7){stop("The threshold number of molecules must be >= 10^7" )}
	#if the above conditions are verified, then the conversion to #cells is done
	#diploid case
	#if(dip){ ncells<-round(quant/6)}
	#haploid case
	#else{ncells<-round(quant/3)}
	#cheking the probabilities input parameters:
	#probEx: extrcation efficiency
	if(!is.numeric(probEx) || is.na(probEx) || probEx <0 || probEx >1)
	{
        stop("'probEx' is a probability, it must belong to [0,1]")
    }
	#probAlq: probability of surviving for aliquots
	if(!is.numeric(probAlq) || is.na(probAlq) || probAlq <0 || probAlq >1)
	{
        stop("'probAlq' is a probability, it must belong to [0,1]")
    }
	
	#probPCR: PCR efficiency
	if(!is.numeric(probPCR) || is.na(probPCR) || probPCR <0 || probPCR >1)
	{
        stop("'pprobPCR' is a probability, it must belong to [0,1]")
    }
	
	#cyc: PCR cycle
	if(!is.numeric(cyc) || is.na(cyc) || cyc <=0)
	{
		stop("'cyc' is the number of PCR cycles, it must be an integer > 0")
	}
	
	#At this point, we have all the input parameters
	#Haploid cells: before the extraction process, we must define the numbers of alleles of type A and B
	
	if(!dip)
	{
		#nAsperm: number of allele of type A, first simulate nAsperm
		nAsperm<-rbinom(n=1,size=ncells,prob=probSperm)
		nBsperm<-ncells-nAsperm
	}
	
	##################1st EXTRACTION STEP
	#dipoloid celles case
	if(dip)
	{
		#alleles of type A surviving the extraction process: nAs, are generated from a binomial distribution
		#with parameters ncells (number of cells) and Probex (extrcation efficiency)
		nAs<- rbinom(n=1,size=ncells,prob=probEx)
		#alleles of type B surviving the extraction process: nBs
		#surviving alleles of type B: nBs
		nBs<- rbinom(n=1,size=ncells,prob=probEx)
	}
	#haploid cells
	else
	{
		nAs<-rbinom(n=1, size=nAsperm,prob=probEx)
		nBs<-rbinom(n=1,size=nBsperm,prob=probEx)
	}
	##################2nd EXTRACTION STEP: aliquots
	#aliquots of type A
	nA<-rbinom(n=1,size=nAs,prob=probAlq)
	#aliquots of type B
	nB<-rbinom(n=1,size=nBs,prob=probAlq)
	
	
	##################PCR efficiency: 
	#for each cycle (defind in cyc), and each allele type
	tmpA<-nA
	tmpB<-nB
	for(t in 1:cyc)
	{
		tmpA<-tmpA + rbinom(n=1,tmpA,prob=probPCR)
		tmpB<-tmpB + rbinom(n=1,tmpB,prob=probPCR)
	}
	#detection threshold T=2x10^7
	#converting from number of molecules to peak heights: to be improved
	#diploid case
	#generating peak heights: this might be subject to change during model calibration
	vecH1 <- round(log((tmpA+Tdrop)/Tdrop) * KH)
    vecH2 <- round(log((tmpB+Tdrop)/Tdrop) * KH)
    a1<-as.integer(vecH1>50)
    a2<-as.integer(vecH2>50)
	n1<-as.integer(tmpA > Tdrop)*tmpA
	n2<-as.integer(tmpB>Tdrop)*tmpB
	Hb<-min(n1,n2)/max(n1,n2)
	if(all(is.na(Hb))){
	Hb<-0
	}
    res <- cbind.data.frame(vecH1, a1, vecH2, a2,Hb)
    colnames(res) <- c("HeightA", "DropA", "HeightB", "DropB",'Hb')
    res
}
