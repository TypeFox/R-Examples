TranSAM <- function(x,samples1,samples2,minChange=0.2,var_filter=0.01,maxFDR=1,changeStep=0.1,scoreFun="magChange"){

	magChange <- function(x,y){
		(median(x)^2)-(median(y)^2)
	}
	
	dstat <- function(x,y,s0){
		numerator <- mean(x)-mean(y)
		a <- ((1/length(x))+(1/length(y)))/(length(x)+length(y)-2)
		denominator <- sqrt(a*((sd(x)^2)+sd(y)^2)) + s0
		numerator/denominator
	}	

	ngenes <- dim(x)[1]

	filtered_genes <- c(1:ngenes)[apply(x,MARGIN=1,sd)>var_filter]
	cat(paste(length(filtered_genes),"genes pass minimum-variance filter \n"))

	if(scoreFun=="dstat"){
		s0 <- as.numeric(quantile(c(apply(x[filtered_genes,samples1],MARGIN=1,sd),apply(x[filtered_genes,samples2],MARGIN=1,sd)),probs=0.05))
		cat(paste("applying SAM d-statistic with fudge-factor",s0," \n"))
	}

	listFun <- function(gene,x,s0,fun,samples1,samples2){
		if(fun=="dstat"){
			score <- dstat(x=x[gene,samples1],y=x[gene,samples2],s0=s0)
		}
		else{
			score <- magChange(x=x[gene,samples1],y=x[gene,samples2])
		}
		score
	}

	realchanges <- unlist(lapply(filtered_genes,listFun,x=x,s0=s0,fun=scoreFun,samples1=samples1,samples2=samples2))
	neg <- which(realchanges < 0)
	realchanges <- sqrt(abs(realchanges))
	realchanges[neg] <- -realchanges[neg]	
	cat(paste("replacing",sum(is.na(realchanges)),"NA values with 0s \n"))
	realchanges[is.na(realchanges)] <- 0

	cat("quantiles of magnitude-changes: \n")
	cat(quantile(realchanges))
	cat("\n")

	realchange_order <- sort(abs(realchanges),decreasing=TRUE,index.return=TRUE)$ix

	# construct balanced permutations of the sample-groups
	cat("creating balanced permutations of sample groups \n")
	nPermSamples <- (floor(min(c(length(samples1),length(samples2)))/2))
	
	combArray1 <- combinations(length(samples1),nPermSamples,v=samples1)
	combArray2 <- combinations(length(samples2),nPermSamples,v=samples2)

	sampleArray <- array(dim=c(dim(combArray1)[1]*dim(combArray2)[1],nPermSamples*2))
	for(i in 1:dim(combArray1)[1]){
		for(j in 1:dim(combArray2)[1]){
			sampleArray[((i-1)*dim(combArray2)[1])+j,] <- c(combArray1[i,],combArray2[j,])
		}	
	}

	comparisonArray <- c()
	for(i in 1:(nrow(sampleArray)-1)){
		for(j in (i+1):nrow(sampleArray)){
			if(sum(sampleArray[i,] %in% sampleArray[j,])==0){
				comparisonArray <- rbind(comparisonArray,c(sampleArray[i,],sampleArray[j,]))
			}
		}
	}

	sampleArray <- unique(comparisonArray)

	groupSel <- c(rep(1,nPermSamples*2),rep(2,nPermSamples*2))
	permutationGroups <- list()
	for(i in 1:dim(sampleArray)[1]){
		permutationGroups[[i]] <- list(grp1=sampleArray[i,groupSel==1],grp2=sampleArray[i,groupSel==2])
	}

	# now calculate changes for every balanced permutation-group
	
	cat(paste("calculating magnitude-changes for",length(permutationGroups),"permutations \n"))
	perm_changes <- array(dim=c(length(filtered_genes),length(permutationGroups)))
	for(i in 1:length(permutationGroups)){
		newchanges <- unlist(lapply(filtered_genes,listFun,x=x,s0=s0,fun=scoreFun,samples1=permutationGroups[[i]]$grp1,samples2=permutationGroups[[i]]$grp2))
		neg <- which(newchanges < 0)
		newchanges <- sqrt(abs(newchanges))
		newchanges[neg] <- -newchanges[neg]
		perm_changes[,i] <- newchanges[realchange_order]
		cat(paste("replacing",sum(is.na(newchanges)),"NA values with 0s \n"))
		newchanges[is.na(newchanges)] <- 0
		# as in SAM, rank the changes in all permutation groups
		if(scoreFun=="dstat"){
			perm_changes[,i] <- sort(newchanges,decreasing=TRUE)
		}
		cat(paste("distribution of permutation",i,"changes: \n"))
		cat(quantile(perm_changes[,i]))
		cat("\n")
	}

	# use these changes to calculate 'expected change' for each gene

	cat("calculating significant genes \n")
	expected_changes <- apply(perm_changes,MARGIN=1,mean)
	cat(paste("distribution of expected changes: \n"))
	cat(quantile(expected_changes))
	cat("\n")
	
	observed_over_expected <- abs(realchanges[realchange_order])-abs(expected_changes)

	FDR <- 1.1
	cat("iterating over thresholds to achieve desired FDR... \n")
	while(FDR>maxFDR){
		significant_genes <- filtered_genes[realchange_order][abs(observed_over_expected) > (minChange)]

		# use table of changes to calculate expected FDR for given threshold-over-expectation
	
		cat(paste("estimating FDR for observed/expected threshold",minChange,"\n"))
		FP_counts <- apply(perm_changes,MARGIN=2,function(x){sum((abs(x)-abs(expected_changes)) > (minChange))})
		FDR <- mean(FP_counts/length(filtered_genes))
	
		cat(paste(length(significant_genes),"genes with changes determined significant by threshold \n"))
		cat(paste("estimated false-discovery rate = ",round(FDR*100,digits=1),"% \n",sep=""))
	
		minChange <- minChange + changeStep
		if(length(significant_genes)==0) stop(paste("no genes found differentially-expressed at FWER",maxFDR))
	}
	
	# output table of genes with most change-over-expected, along with the predicted FDR
	significant_ratios <- observed_over_expected[filtered_genes[realchange_order] %in% significant_genes]
        significant_changes <- realchanges[realchange_order][filtered_genes[realchange_order] %in% significant_genes]
        siggene_ordering <- order(abs(significant_ratios),decreasing=TRUE)
	significant_genes <- significant_genes[siggene_ordering]
	significant_ratios <- significant_ratios[siggene_ordering]
	significant_changes <- significant_changes[siggene_ordering]
	

	cat(paste("estimating q-values for",length(significant_genes),"genes \n"))
	FDRestimates <- rep(NA,length(significant_genes))
	for(i in 1:length(significant_genes)){
		thisChange <- significant_ratios[i]
		this.genes <- filtered_genes[realchange_order][abs(observed_over_expected) >= (thisChange)]
		this.FPcounts <- apply(perm_changes,MARGIN=2,function(x){sum((abs(x)-abs(expected_changes)) >= (thisChange))})
		FDRestimates[i] <- mean(this.FPcounts/length(filtered_genes))
	}

	gene.names <- significant_genes
	if(!is.null(rownames(x))) gene.names <- rownames(x)[significant_genes]	
	cat("outputting table \n")
	data.frame(genes=gene.names,obs.exp.ratios=significant_ratios,change=significant_changes,FDR.estimate=FDRestimates)

}
