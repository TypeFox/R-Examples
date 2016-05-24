## split correlation matrix into two blocks
## centers are points with lowest value
## calculate memberships but (!) with same number of data points per cluster
## therefore split split into halves according to order of correlation
blxsplit <- function(x, bs){
    n <- nrow(x)
    n.nna <- sqrt(sum(!is.na(x)))
    centers <- which(x==min(x, na.rm=T), arr.ind=T)[1,]
    assocs <- order(apply(x[centers,], 2, function(u) ifelse(all(is.na(u)), NA, u[which.max(u)]*sign(which.max(u)-1.5))))[1:n.nna]
    list(assocs[1:(trunc(ceiling(n.nna/bs)/2)*bs)], assocs[(trunc(ceiling(n.nna/bs)/2)*bs+1):n.nna])
}

## recursive split until all blocks are formed
## no tail recursion in R -> ugly code
bloxs <- function(x, bs=4){
    n <- nrow(x)
    if(n <= bs)
        stop("Block size must be smaller than input dimensions.")
    blx.tosplit <- list(1:n)
    blx <- list()
    while(length(blx.tosplit)>0){        
        x.cur <- matrix(NA, nrow=n, ncol=n)
        x.cur[blx.tosplit[[1]], blx.tosplit[[1]]] <- x[blx.tosplit[[1]], blx.tosplit[[1]]]
        blx.cur <- blxsplit(x.cur, bs)
        blx <- c(blx, blx.cur[sapply(blx.cur, length)<=bs])
        blx.tosplit <- c(blx.tosplit[-1], blx.cur[sapply(blx.cur, length)>bs])
    }
    blx
}

evaluateGeneSetUncertainty <- function(
		...,
		dat,
		geneSet,
		analysis,
		numSamplesUncertainty,
		blockSize				= 1,#round(length(geneSet)^(0.25))+1,
		k						= seq(0.01, 0.99, by=0.01),
		signLevel 				= 0.05,
		preprocessGeneSet		= FALSE,
		cluster 				= NULL){

	if(analysis$significance != "significance.sampling"){
		stop("'evaluateGeneSetUncertainty': Use an analysis with significance assessment 'significance.sampling'.")
	}

	##########################################
	# init seeds
	##################
	## adapted from
	## package:		parallel
	## function:	clusterSetRNGStream
	##########################################
	if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
		oldseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	}else{
		oldseed <- NULL
	}
	iseed <- round(runif(n=1) * .Machine$integer.max)
	RNGkind("L'Ecuyer-CMRG")
	set.seed(iseed)

	seedList <- vector(mode = "list", length= length(k))
	seedList[[1]] <- .Random.seed
	if(length(k) > 1){
		for (i in 2:length(k)) {
			seedList[[i]] <- nextRNGStream(seedList[[i-1]])
		}
	}
	##########################################

	if(preprocessGeneSet){
		geneSet <- preprocessGs(dat, list(originalGeneSet = geneSet))
		rownames(dat) <- tolower(rownames(dat))
	}else{
		geneSet <- list(originalGeneSet = geneSet)
	}

	##########################################
	# original gene set
	##########################################
	originalGs <- geneSetAnalysis(
		...,
		dat						= dat,
		geneSets				= geneSet,
		analysis				= analysis,
		signLevel				= signLevel,
		adjustmentMethod		= "none")

	transformation <- originalGs$res.all[[1]]$geneSetValues$transformation
	gls <- originalGs$res.all[[1]]$geneSetValues$gls

	if(blockSize > 1){	
		datGeneSet <- dat[rownames(dat)%in%unlist(geneSet),]
		geneSetCorrelations <- cor(t(datGeneSet), t(datGeneSet))
		blockOrder <- unlist(bloxs(geneSetCorrelations, blockSize))
	}else{
		blockOrder <- 1:length(unlist(geneSet))
	}
	geneSet <- unlist(geneSet)[blockOrder]

	##########################################
	# parallel: yes/no
	##########################################
	if(is.null(cluster)){
		myApply <- mapply
	}else{
		clusterCall(cluster, function() library("GiANT", character.only=TRUE))
		clusterExport(cluster,
			c("dat","geneSet","analysis","args","signLevel"),
			envir = environment())
		#export functions possibly defined by the user
		clusterExport(cluster, c(analysis$gls,
			analysis$transformation,
			analysis$gss,
			analysis$globalStat,
			analysis$significance))

		myApply <- function(...){
			clusterMap(cl = cluster, ...)
		}
	}

	##########################################
	# values for #numSamplesUncertainty different
	# (partially) random genesets
	##########################################
	if(!is.null(cluster)){
		clusterExport(cluster,
			c("gls", "transformation", "originalGs"),
			envir = environment())
	}

	gs.stat <- myApply(function(x,s){
			assign(".Random.seed", s, envir = .GlobalEnv)
			
			# get names which are not in the gene set
			selectable <- rownames(dat)[!rownames(dat)%in%geneSet]

			if(blockSize >= length(geneSet)){
				stop("blockSize (", blockSize,") >= geneSet size (", length(geneSet), ")")
			}

			genesToSample <- round((1-x)*length(geneSet))

			genesetBlocks <- seq(1,length(geneSet),blockSize)
			genesetBlockSizes <- c(genesetBlocks[2:length(genesetBlocks)], length(geneSet)+1)-genesetBlocks

			allgenesets <- replicate(numSamplesUncertainty, {

					blocksToReplace <- sample(length(genesetBlocks), length(genesetBlocks))
					genesToReplace <- unlist(lapply(blocksToReplace, function(i){
							return(genesetBlocks[i]:(genesetBlocks[i]+min(c((genesetBlockSizes[i]-1), length(geneSet)))))
						}))
					genesToReplace <- genesToReplace[1:genesToSample] #NEW
					replacingGenes <- sample(selectable, genesToSample)#replacingGenes[1:length(genesToReplace)]
					newGeneset <- geneSet
					newGeneset[genesToReplace] <- replacingGenes #selectable[replacingGenes]

					return(newGeneset)
				})
			
			colnames(allgenesets) <- paste("uncertainGeneSet",1:numSamplesUncertainty,sep ="")
			rownames(allgenesets) <- NULL

			gssValues <- apply(allgenesets, 2, doGSS,
				dat = dat,
				analysis = analysis,
				parameters = list(...),
				transformation = transformation)

			names(gssValues) <- paste("iter",1:numSamplesUncertainty,sep ="")

			confidenceInterval <- quantile(c(gssValues,originalGs$res.all[[1]]$geneSetValues$gss), probs = c(signLevel, 0.5, 1-signLevel))
			names(confidenceInterval) <- paste(c(signLevel, 0.5, 1-signLevel)*100, "%-quantile", sep ="")

			return(list(
				confidenceValues = confidenceInterval,
				gssValues = gssValues,
				uncertainGeneSets = allgenesets,
				k = x))
		}, k, seedList, SIMPLIFY=FALSE)

	names(gs.stat) <- paste(k*100, "%-originalGenes", sep ="")

	conf <- t(sapply(gs.stat, "[[", 1))
	rownames(conf) <- paste(k*100, "%-originalGenes", sep ="")

	quant <- c(signLevel, 0.5, 1-signLevel)

	nullDistr <- quantile(originalGs$res.all[[1]]$significanceValues$gssValues, probs = quant)

	ind <- max(which(abs(conf[,1]-nullDistr[3]) == min(abs(conf[,1]-nullDistr[3]))))
	if(originalGs$analysis$testAlternative == "greater"){
		if((conf[,1]-nullDistr[3])[ind] > 0){
			uncertainty <- k[ind]
		}else{
			uncertainty <- k[ind+1]
		}
		#uncertainty <- k[which((abs(conf[,1]-nullDistr[3]) == min(abs(conf[,1]-nullDistr[3]))) & (conf[,1]-nullDistr[3]) > 0)]
	}else if(originalGs$analysis$testAlternative == "less"){
		if((conf[,1]-nullDistr[3])[ind] < 0){
			uncertainty <- 1-k[ind]
		}else{
			uncertainty <- 1-k[ind+1]
		}
		#uncertainty <- k[which((abs(conf[,3]-nullDistr[1]) == min(abs(conf[,3]-nullDistr[1]))) & (conf[,3]-nullDistr[1]) < 0)]
	}

	sstat <- list(
		uncertainty = uncertainty,
		confidenceValues = conf,
		uncertaintyEvaluations = gs.stat,
		signLevel = signLevel,
		originalGeneSetValues = originalGs)

	class(sstat) <- "uncertaintyResult"
	return(sstat)
}
