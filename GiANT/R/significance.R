significance.sampling <- function(
		...,
		dat,
		geneSet,
		analysis,
		glsValues,
		numSamples = 1000){

	#################################################
	#sampele random gene sets
	randomGeneSets <- replicate(numSamples, sample(x=rownames(dat), size=length(geneSet)))

	#gss for all randomly sampeled gene sets
	gssValues <- apply(randomGeneSets, 2, doGSS,
		dat = dat,
		analysis = analysis,
		parameters = list(...),
		transformation = glsValues)
	names(gssValues) <- paste("iter",1:numSamples,sep="")

	sstat <- list(gssValues = gssValues, randomGeneSets = randomGeneSets)
	return(sstat)
}

#####
# Not required when all genes (or all blocks) in a gene set are
# replaced with random genes.
# Can be used when also random genes are sampeled as blocks.
# Currently a computation of pairwise correlations for
# more than 10000 genes is not practical.
#
# significance.blockResampling <- function(
# 		...,
# 		dat,
# 		geneSet,
# 		analysis,
# 		glsValues,
# 		numSamples = 1000,
# 		blockSize = 1){#round(length(geneSet)^(0.25))+1

# 	#################################################
# 	#sample random gene sets
# 	#randomGeneSets <- replicate(numSamples, sample(x=rownames(dat), size=length(geneSet)))

# 	nms <- rownames(dat)[order(glsValues)]

# 	selectable <- nms[!nms%in%geneSet]

# 	genesetBlocks <- seq(1,length(geneSet),blockSize)
# 	genesetBlockSizes <- c(genesetBlocks[2:length(genesetBlocks)], length(geneSet)+1)-genesetBlocks

# 	randomGeneSets <- replicate(numSamples, {
# 			blocksToReplace <- sample(length(genesetBlocks), length(genesetBlocks))

# 			genesToReplace <- unlist(lapply(blocksToReplace, function(i){
# 					return(genesetBlocks[i]:(genesetBlocks[i]+min(c((genesetBlockSizes[i]-1), length(geneSet)))))
# 				}))

# 			genesToReplace <- genesToReplace[1:length(geneSet)]

# 			replacingGenes <- sample(selectable, length(geneSet))

# 			newGeneset <- geneSet
# 			newGeneset[genesToReplace] <- replacingGenes

# 			return(newGeneset)
# 		})

# 	#gss for all randomly sampeled gene sets
# 	gssValues <- apply(randomGeneSets, 2, doGSS,
# 		dat = dat,
# 		analysis = analysis,
# 		parameters = list(...),
# 		transformation = glsValues)
# 	names(gssValues) <- paste("iter",1:numSamples,sep="")

# 	sstat <- list(gssValues = gssValues, randomGeneSets = randomGeneSets)
# 	return(sstat)
# }

significance.permutation <- function(
		...,
		dat,
		geneSet,
		analysis,
		glsValues,
		numSamples = 1000,
		labs){

	parameters <- list(...)

	cl <- unique(labs)
	basic <- rep(cl[2], length(labs))

	geneSetIndices <- which(tolower(rownames(dat)) %in% tolower(geneSet))

	#Calculate the number of possible permutations for the given number of samples
	#and check whether this number is bigger than numSamples or not.
	nPossiblePermutations <- choose(length(labs), sum(labs == cl[1]))

	if(nPossiblePermutations < numSamples){
		warning(paste("Number of possible permutations (",
			nPossiblePermutations,
			") is smaller than given number of numSamples (",
			numSamples,
			").\nCalculating all possible permutations for testing.", sep=""))
		
		permutations <- combn(length(labs), sum(labs == cl[1]))
		numSamples <- nPossiblePermutations
	}else{
		permutations <- replicate(numSamples, sample(length(labs),sum(labs == cl[1])))
	}
	colnames(permutations) <- paste("iter",1:numSamples,sep="")
	
	#gene level statistic
	originalParameter <- labs

	glsValues <- apply(permutations, 2, function(x){
			rndPerm <- basic
			rndPerm[x] <- cl[1]

			parameters[["labs"]] <- rndPerm

			doGLS(dat = dat,
				analysis = analysis,
				parameters = parameters)
		})
	rownames(glsValues) <- rownames(dat)
	colnames(glsValues) <- paste("iter",1:numSamples,sep="")

	parameters[["labs"]] <- originalParameter

	#transformation
	transformedValues <- apply(glsValues, 2, function(x){
			doTransformation(
				analysis = analysis,
				parameters = parameters,
				gls = x)
		})
	rownames(transformedValues) <- rownames(dat)
	colnames(transformedValues) <- paste("iter",1:numSamples,sep="")

	#gene set statistic
	gssValues <- apply(transformedValues, 2, function(x){
			doGSS(dat = dat,
				geneSet = geneSet,
				analysis = analysis,
				parameters = parameters,
				transformation = x)
		})
	names(gssValues) <- paste("iter",1:numSamples,sep="")

	sstat <- list(gssValues = gssValues, permutations = permutations)
	return(sstat)
}

significance.restandardization <- function(
		...,
		dat,
		geneSet,
		analysis,
		glsValues,
		numSamples = 1000,
		labs){
	#resampling of genes
	samp <- significance.sampling(
		...,
		dat = dat,
		geneSet = geneSet,
		analysis = analysis,
		numSamples = numSamples)

	mean_samp <- mean(samp$gssValues)
	sd_samp <- sd(samp$gssValues)

	#permutation of labels
	perm <- significance.permutation(
		...,
		dat = dat,
		geneSet = geneSet,
		analysis = analysis,
		numSamples = numSamples,
		labs = labs)

	mean_perm <- mean(perm$gssValues)
	sd_perm <- sd(perm$gssValues)

	#combine values of permutation test statistic with resampling values
	resta <- mean_samp + ((sd_samp/sd_perm)*((perm$gssValues - mean_perm)/sqrt(length(geneSet))))

	sstat <- list(gssValues = resta, samplingValues = samp, permutationValues = perm)
}