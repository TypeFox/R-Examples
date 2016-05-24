burdenMeta <- function(..., SNPInfo=NULL, wts = 1, snpNames = "Name", aggregateBy = "gene", mafRange = c(0,0.5), verbose=FALSE){
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	}	
	
	cohorts <- list(...)	
	ncohort <- length(cohorts)
	if(any(unlist(lapply(cohorts,class)) != "skatCohort")){
	 	stop("argument to ... is not a skatCohort object!")
	}
	
	genelist <- na.omit(unique(SNPInfo[,aggregateBy]))
	res.strings <- data.frame("gene"=genelist,stringsAsFactors=F)
	res.numeric <- matrix(NA, nrow= nrow(res.strings),ncol =  length(c("p","beta","se","cmafTotal", "cmafUsed","nsnpsTotal","nsnpsUsed","nmiss")))
	colnames(res.numeric) <- c("p","beta","se","cmafTotal", "cmafUsed","nsnpsTotal","nsnpsUsed","nmiss")

	if(verbose){
    	cat("\n Meta Analyzing... Progress:\n")
    	pb <- txtProgressBar(min = 0, max = length(genelist), style = 3)
    	pb.i <- 0
    }

	snp.names.list <- split(SNPInfo[,snpNames],SNPInfo[,aggregateBy])
	
	ri <- 0
	for(gene in genelist){
		
		ri <- ri+1
		nsnps.sub <- length(snp.names.list[[gene]])
		
		maf <- vcount <- numeric(nsnps.sub)
		n.total <- numeric(nsnps.sub)
		n.miss <- numeric(nsnps.sub)
		big.cov <- matrix(0, nsnps.sub,nsnps.sub)
		mscores <- numeric(nsnps.sub)
		
		vary.ave <- 0
		for(cohort.k in 1:ncohort){
			cohort.gene <- cohorts[[cohort.k]][[gene]]
			
			if(!is.null(cohort.gene)){
				cohort.gene <- lapply(cohort.gene,function(x){replace(x,is.nan(x),0)})
				sub <- match(snp.names.list[[gene]],colnames(cohort.gene$cov))
				if(any(is.na(sub)) | any(sub != 1:length(sub), na.rm=TRUE) | length(cohort.gene$maf) > nsnps.sub){
						cohort.gene$maf <- cohort.gene$maf[sub]
						cohort.gene$maf[is.na(sub)] <- -1
						
						cohort.gene$cov <- cohort.gene$cov[sub,sub,drop=FALSE]
						cohort.gene$cov[is.na(sub),] <- cohort.gene$cov[,is.na(sub)] <- 0
												
						cohort.gene$scores <- cohort.gene$scores[sub]
						cohort.gene$scores[is.na(sub)] <- 0
				}
				#imputation
				n.miss[cohort.gene$maf < 0] <- n.miss[cohort.gene$maf < 0] + cohort.gene$n
				n.total[cohort.gene$maf >= 0] <- n.total[cohort.gene$maf >= 0]+cohort.gene$n
				cohort.gene$maf[cohort.gene$maf < 0] <- 0
				
				maf <- maf + 2*cohort.gene$maf*cohort.gene$n
				mscores <- mscores + cohort.gene$scores/cohort.gene$sey^2
				big.cov <- big.cov + cohort.gene$cov/cohort.gene$sey^2
				vary.ave <- vary.ave + max(cohort.gene$n,na.rm=T)*cohort.gene$sey^2
			}
		}

		maf <- maf/(2*n.total)
		maf[is.nan(maf)] <- 0
		maf <- sapply(maf, function(x){min(x,1-x)})
		
		if(is.function(wts)){
			tmpwts <- ifelse(maf > 0, wts(maf),0)
		} else if(is.character(wts)){
			tmpwts <- as.numeric(SNPInfo[SNPInfo[,aggregateBy]==gene,wts])
		} else {
			tmpwts <- rep(1,length(maf))
		}
		
		if( !all(mafRange == c(0,0.5))){
			keep <- (maf > min(mafRange)) & (maf <= max(mafRange))
			tmpwts[!keep] <- 0
		}
		
		vary.ave <- vary.ave/max(n.total,na.rm=TRUE)
		bscorevar <- c(tmpwts%*%big.cov%*%tmpwts)
		bscore <- tmpwts%*%mscores 
		
		res.numeric[ri,"beta"] = ifelse(bscorevar !=0, bscore/bscorevar, NA)
		res.numeric[ri,"se"] = sqrt(1/bscorevar)
		res.numeric[ri,"cmafTotal"] = sum(maf,na.rm=TRUE)
		res.numeric[ri,"cmafUsed"] = sum(maf[tmpwts != 0],na.rm=TRUE)
		res.numeric[ri,"nsnpsTotal"] = length(maf)
		res.numeric[ri,"nmiss"] = sum(n.miss[tmpwts != 0], na.rm =T)
		res.numeric[ri,"nsnpsUsed"] = sum(tmpwts != 0)
		res.numeric[ri,"p"] = ifelse(bscorevar !=0,pchisq(bscore^2/bscorevar,lower.tail=FALSE,df=1),NA)
		if(verbose){
			pb.i <- pb.i+1
			setTxtProgressBar(pb, pb.i)
		}
	
	}
	if(verbose) close(pb)
	return(cbind(res.strings,res.numeric))
}
