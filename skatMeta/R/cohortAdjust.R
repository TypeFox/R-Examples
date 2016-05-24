skatCohortAdj <- function(Z, formula, family = gaussian(), SNPInfo=NULL, adjustments= NULL, snpNames = "Name", aggregateBy = "gene", data=parent.frame()){
	if(is.null(SNPInfo)) stop("SNPInfo file must be provided!")
	if(is.null(adjustments)) stop("adjustments must be provided!")
	
	re <- NULL
	#update formula
	null.formula.adj <- update(formula, ".~.+ geno.adj")
	
	genes <- unique(adjustments[,aggregateBy])
	for(gene in genes){
		#subset to the gene:
		minisnpinfo <- SNPInfo[SNPInfo[,aggregateBy] == gene,]
		
		#remove adjustment snp
		adjustmentsnps <- subset(adjustments[,snpNames],adjustments[,aggregateBy] ==gene)
		minisnpinfo <- subset(minisnpinfo, !(minisnpinfo[,snpNames] %in% adjustmentsnps) , drop=FALSE)
		
		#if there are more snps in this gene....
		if(nrow(minisnpinfo) > 0){
			#add genotype to adjustment variables
			zz <- subset(Z,select = colnames(Z) %in% adjustmentsnps,drop=FALSE)
			zz <- apply(zz,2,function(z){
				z[is.na(z)] <- mean(z,na.rm=T)
				z
			})
			data$geno.adj <- zz
			
			Z.sub <- subset(Z, select = colnames(Z) %in% minisnpinfo$Name,drop=FALSE)
			
			#if our sample has more snps in this gene...
			if(ncol(Z.sub) > 0){	
				re.tmp <- skatCohort(Z= Z.sub, formula = null.formula.adj, family = family, 
					SNPInfo = minisnpinfo, aggregateBy = aggregateBy, snpNames = snpNames, data=data)
				class(re.tmp) <- "list"
				if(is.null(re)){
					re <- re.tmp
				} else {
					re <- c(re,re.tmp)
				}		
			} 
		}	
	}
	class(re) <- "skatCohort"
	return(re)
}


skatFamCohortAdj <- function(Z, formula, SNPInfo=NULL, adjustments= NULL, snpNames = "Name", aggregateBy = "gene",fullkins, sparse = TRUE, data=parent.frame()){
	if(is.null(SNPInfo)) stop("SNPInfo file must be provided!")
	if(is.null(adjustments)) stop("adjustments must be provided!")
	
	re <- NULL
	#update formula
	null.formula.adj <- update(formula, ".~.+ geno.adj")
	
	genes <- unique(adjustments[,aggregateBy])
	for(gene in genes){
		#subset to the gene:
		minisnpinfo <- SNPInfo[SNPInfo[,aggregateBy] == gene,]
		
		#remove adjustment snp
		adjustmentsnps <- subset(adjustments[,snpNames],adjustments[,aggregateBy] ==gene)
		minisnpinfo <- subset(minisnpinfo, !(minisnpinfo[,snpNames] %in% adjustmentsnps) , drop=FALSE)
		
		#if there are more snps in this gene....
		if(nrow(minisnpinfo) > 0){
			#add genotype to adjustment variables
			zz <- subset(Z,select = colnames(Z) %in% adjustmentsnps,drop=FALSE)
			zz <- apply(zz,2,function(z){
				z[is.na(z)] <- mean(z,na.rm=T)
				z
			})
			data$geno.adj <- zz
			
			Z.sub <- subset(Z, select = colnames(Z) %in% minisnpinfo$Name,drop=FALSE)			
			#if our sample has more snps in this gene...
			if(ncol(Z.sub) > 0){	
				re.tmp <- skatFamCohort(Z= Z.sub, formula = null.formula.adj, 
					SNPInfo = minisnpinfo, aggregateBy = aggregateBy, snpNames = snpNames, 
						fullkins, sparse = sparse, data=data)
				class(re.tmp) <- "list"
				if(is.null(re)){
					re <- re.tmp
				} else {
					re <- c(re,re.tmp)
				}		
			} 
		}	
	}
	class(re) <- "skatCohort"
	return(re)
}