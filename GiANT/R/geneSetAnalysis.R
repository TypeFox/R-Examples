##########################################################################
#geneSetAnalysis
#
# dat:				A matrix of genes expressions from different samples
#					(rows = genes).
# geneSets:			A list of gene sets (vectors of gene names) each gene
#					should be present in rownames(dat).
# analysis:			The analysis which should be performed, a gsAnalysis
#					object.
# signLevel:		significance level for the applied tests.
# preprocessGeneSets:Preprocess gene sets
#						- set all entries to lower case
#						- remove genes not in rownames(dat)
# adjustmentMethod:	The method used by p.adjust to correct for
#					multiple testing (one method from p.adjust.methods).
# cluster:			If an initialised cluster is supplied gene sets are
#					analysed in parallel.
##########################################################################
geneSetAnalysis <- function(
		...,
		dat,
		geneSets,
		analysis,
		signLevel 				= 0.05,
		preprocessGeneSets		= FALSE,
		adjustmentMethod		= p.adjust.methods,
		cluster 				= NULL){

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

	seedList <- vector(mode = "list", length= length(geneSets))
	seedList[[1]] <- .Random.seed
	if(length(geneSets) > 1){
		for (i in 2:length(geneSets)) {
			seedList[[i]] <- nextRNGStream(seedList[[i-1]])
		}
	}
	##########################################

	adjustmentMethod <- match.arg(adjustmentMethod)

	if(class(analysis)!="gsAnalysis"){
		stop("Analysis must be of class gsAnalysis.\n
			Use function 'gsAnalysis' to generate an object of this class.")
	}

	if(preprocessGeneSets){
		geneSets <- preprocessGs(dat, geneSets)
		rownames(dat) <- tolower(rownames(dat))
	}

	gs_names <- names(geneSets)

	if(is.null(gs_names)){
		gs_names <- paste("geneset",1:length(geneSets),sep="")
		names(geneSets) <- gs_names
	}

	args <- list(...)
	analysisType <- ""

	##########################################
	# parallel: yes/no
	##########################################
	if(is.null(cluster)){
		myApply <- mapply
	}else{
		clusterCall(cluster, function() library("GiANT", character.only=TRUE))
		clusterExport(cluster,
			c("dat","geneSets","analysis","args","signLevel"),
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
	# analyses for each geneset
	##########################################
	if(!is.null(analysis$globalStat)){
		##########################################
		# global analyses
		##########################################
		gs.stat <- myApply(function(x, s){
				assign(".Random.seed", s, envir = .GlobalEnv)
				#################################################
				#apply global analysis
				gs <- doGlobalAnalysis(
					dat = dat,
					geneSet = x,
					analysis = analysis,
					parameters = args)
				#################################################
				return(list(pValue = gs$pValue,
					geneSetValues = gs,
					geneSet = x))
			}, geneSets, seedList, SIMPLIFY=FALSE)
		analysisType <- "globalAnalysis"
	}else{
		##########################################
		# enrichment analysis procedure
		# gls -> transformation -> gss
		##########################################
		#gene level statistic
		gls <- doGLS(
			dat = dat,
			analysis = analysis,
			parameters = args)
		##########################################
		#transformation
		transformation <- doTransformation(
			analysis = analysis,
			parameters = args,
			gls = gls)

		if(!is.null(cluster)){
			clusterExport(cluster,
				c("gls", "transformation"),
				envir = environment())
		}

		gs.stat <- myApply(function(x,s){
				assign(".Random.seed", s, envir = .GlobalEnv)
				#################################################
				#gene set statistic
				gss <- doGSS(
					dat = dat,
					geneSet = x,
					analysis = analysis,
					parameters = args,
					transformation = transformation)
				#################################################
				#significance assessment
				gsSig <- doSignificanceAssessment(
					dat = dat,
					geneSet = x,
					analysis = analysis,
					glsValues = transformation,
					parameters = args)
				#################################################
				# calc p-Value
				if(is.list(gss) && is.null(gsSig)){
					return(list(pValue = gss$pValue,
						geneSetValues = c(gss, list(gls = gls, transformation = transformation)),
						geneSet = x))
				}else if(is.null(gsSig)){
					pValue <- NULL
				}else{
					if(analysis$testAlternative == "greater"){
						pValue <- sum(c(gsSig$gssValues, gss) >= gss) / length(c(gsSig$gssValues, gss))
					}else if(analysis$testAlternative == "less"){
						pValue <- sum(c(gsSig$gssValues, gss) <= gss) / length(c(gsSig$gssValues, gss))
					}
					# ??? two.sided ???
					#else{
					#	pValue <- mean(c(sum(c(gsSig$gssValues, gss) >= gss) / length(c(gsSig$gssValues, gss)),
					#		sum(c(gsSig$gssValues, gss) <= gss) / length(c(gsSig$gssValues, gss))))
					#}
				}
				#################################################
				return(list(pValue = pValue,
					geneSetValues = list(gls = gls, transformation = transformation, gss = gss),
					significanceValues = gsSig,
					geneSet = x))
			}, geneSets, seedList, SIMPLIFY=FALSE)
		analysisType <- "geneSetAnalysis"
	}

	names(gs.stat) <- gs_names
	#################################################
	# correction for multiple testing
	rawPValues <- sapply(gs.stat, "[[", 1)
	adjustedPValues <- p.adjust(rawPValues, method = adjustmentMethod)
	#################################################
	if(is.null(oldseed)){
		rm(.Random.seed, envir = .GlobalEnv)
	}else{
		assign(".Random.seed", oldseed, envir = .GlobalEnv)
	}

	result <- list(
		adjustedPValues = adjustedPValues,
		rawPValues = rawPValues,
		res.all = gs.stat,
		signLevel = signLevel,
		analysis = analysis,
		analysisType = analysisType,
		adjustmentMethod = adjustmentMethod)
	class(result) <- "gsaResult"
	return(result)
}
