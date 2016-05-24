singlesnpMeta <- function(..., SNPInfo=NULL, snpNames = "Name", aggregateBy = "gene", cohortBetas = TRUE, verbose = FALSE){
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	}
	cohorts <- list(...)	
	ncohort <- length(cohorts)
	cl <- match.call(expand.dots=FALSE)
	cohort.names <- NULL
	for(i in 1:ncohort) cohort.names <- c(cohort.names,as.character(cl[[2]][[i]]))
	if(any(unlist(lapply(cohorts,class)) != "skatCohort")){
	 	stop("argument to ... is not a skatCohort object!")
	}
	
	genelist <- na.omit(unique(SNPInfo[,aggregateBy]))
	res.strings <- data.frame("gene"=as.character(SNPInfo[,aggregateBy]),
					"Name" = as.character(SNPInfo[,snpNames]),stringsAsFactors=F)
	res.numeric <- matrix(NA, nrow= nrow(res.strings),ncol =  length(c("p","maf","nmiss","ntotal", "beta" ,"se" )))
	colnames(res.numeric) <- c("p","maf","nmiss","ntotal", "beta" ,"se" )	
					

	if(cohortBetas){
		resdf.cohort = matrix(NA,nrow=nrow(SNPInfo),ncol=2*ncohort)
		oo <- as.numeric(gl(ncohort,2))
		oo[1:ncohort*2]<- oo[1:ncohort*2]+ncohort
		colnames(resdf.cohort) <- c(paste(c("beta"),cohort.names,sep="."),paste(c("se"),cohort.names ,sep="."))[oo]
	}
	
	if(verbose){
    	cat("\n Meta Analyzing... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = length(genelist), style = 3)
    	pb.i <- 0
    }

	#gene.namelist <- split(SNPInfo[,snpNames],SNPInfo[,aggregateBy])
	ris <- split(1:nrow(SNPInfo),SNPInfo[,aggregateBy])
	for(gene in genelist){
		
		ri <- ris[[gene]]
		SNPInfo.sub <- SNPInfo[ri,,drop =FALSE]
		nsnps.sub <- nrow(SNPInfo.sub)
		
		maf <- vcount <- scorevar <- mscore <- numeric(nsnps.sub)
		
		n.total <- numeric(nsnps.sub)
		n.miss <- numeric(nsnps.sub)
		vary.ave <- 0
		for(cohort.k in 1:ncohort){
			cohort.gene <- cohorts[[cohort.k]][[gene]]
			
			if(!is.null(cohort.gene)){
				cohort.gene <- lapply(cohort.gene,function(x){replace(x,is.nan(x),0)})
				sub <- match(SNPInfo.sub[,snpNames],colnames(cohort.gene$cov))
				if(any(is.na(sub)) | any(sub != 1:length(sub), na.rm=TRUE) | length(cohort.gene$maf) > nsnps.sub){
						if(any(is.na(sub))) warning("Some SNPs were not in SNPInfo file for gene ", gene," and cohort ",cohort.names[[cohort.k]])							
						cohort.gene$maf <- cohort.gene$maf[sub]
						cohort.gene$maf[is.na(sub)] <- -1
							
						cohort.gene$scores <- cohort.gene$scores[sub]
						cohort.gene$scores[is.na(sub)] <- 0
							
						cohort.gene$cov <- cohort.gene$cov[sub,sub,drop=FALSE]
						cohort.gene$cov[is.na(sub),] <- cohort.gene$cov[,is.na(sub)] <- 0		
					}				
					
					n.total <- n.total + (cohort.gene$maf >= 0)*cohort.gene$n				
					n.miss[cohort.gene$maf < 0] <- n.miss[cohort.gene$maf < 0] + cohort.gene$n
					cohort.gene$maf[cohort.gene$maf < 0] <- 0
	
					mscore <- mscore + cohort.gene$score/cohort.gene$sey^2	
					scorevar <- scorevar + diag(cohort.gene$cov)/cohort.gene$sey^2
					maf <- maf + 2*cohort.gene$maf*(cohort.gene$n)
					vary.ave <- vary.ave + max(cohort.gene$n,na.rm=T)*cohort.gene$sey^2
		
					if(cohortBetas){
						resdf.cohort[ri,paste(c("beta"),cohort.names[cohort.k] ,sep=".")] <- cohort.gene$score/diag(cohort.gene$cov)
						resdf.cohort[ri,paste(c("se"), cohort.names[cohort.k] ,sep=".")] <- cohort.gene$sey/sqrt(diag(cohort.gene$cov))
					}
				}
		}
		vary.ave <- vary.ave/n.total

		maf <- maf/(2*n.total)
		maf[is.nan(maf)] <- 0
		maf <- sapply(maf, function(x){min(x,1-x)})
				
		res.numeric[ri,c("beta","se","maf","nmiss","ntotal","p")] <- cbind( ifelse(scorevar !=0, mscore/scorevar, NA),
																	sqrt(1/scorevar),
																	maf,
																	n.miss,
																	n.total,
																	ifelse(scorevar !=0, pchisq(mscore^2/scorevar,lower.tail=FALSE,df=1), NA))

		if(verbose){
			pb.i <- pb.i+1
			setTxtProgressBar(pb, pb.i)
		}
	}
	if(verbose) close(pb)
	
	if(cohortBetas){
		return(cbind(res.strings,res.numeric,resdf.cohort))
	} else {
		return(cbind(res.strings,res.numeric))
	}
}
