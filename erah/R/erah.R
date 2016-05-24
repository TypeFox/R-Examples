
globalVariables("mslib")

## MetaboSet Class Definition:

setClass(Class = "eRahSoftParameters", representation = representation(algorithm="character", min.peak.width = "numeric", min.peak.height = "numeric", noise.threshold = "numeric", avoid.processing.mz = "vector",compression.coef = "numeric", analysis.time="vector"))

setClass(Class = "eRahAlParameters", representation = representation(algorithm="character", min.spectra.cor="numeric", max.time.dist="numeric", mz.range="vector", method="character"))	

setClass(Class = "eRahIdParameters", representation = representation(algorithm="character", database="character", compare.only.mz="vector"))

setClass(Class = "MetaData", representation = representation(Instrumental = "data.frame", Phenotype = "data.frame", DataDirectory="character"))

setClass(Class = "Statistics", representation = representation(Univariate="data.frame", Multivariate="data.frame"))	

setClass(Class="MSResultsParameters", representation = representation(Alignment = "list", Identification = "list"))

setClass(Class="Data", representation = representation(FeatureList = "list", FactorList = "list", Parameters = "list"))

setClass(Class = "Results", representation = representation(Parameters="MSResultsParameters", Alignment = "data.frame", Identification="data.frame", Statistics="Statistics"))

setClass(Class="MetaboSet",representation= representation(Info = "character", Data="Data", MetaData="MetaData", Results = "Results"))


## Intern Classes for eRah:

setClass(Class = "eRah_DB", representation = representation(name="character", version="character", info="character", database="list"))

setClass(Class = "RawDataParameters", representation = representation(data = "matrix", min.mz = "numeric", max.mz = "numeric", start.time = "numeric", mz.resolution = "numeric", scans.per.second = "numeric", avoid.processing.mz = "vector", min.peak.width = "numeric", min.peak.height = "numeric", noise.threshold = "numeric", compression.coef = "numeric"))
	
	
setMethod("show", "MetaboSet", function(object){
	cat("A \"MetaboSet\" object containing", length(object@Data@FactorList), "samples \n \n" ,sep=" ")
	cat("Data processed with", object@Data@Parameters$algorithm, "\n" ,sep=" ")
	cat("Info attached for this experiment: \n", object@Info)
})

metaData <- function(object) {object@MetaData@Instrumental}
phenoData <- function(object) {object@MetaData@Phenotype}
		
	
#setClass(Class="expClasses",representation= representation(classes.type = "character", classes.summary = "data.frame"))

# setMethod("show", "expClasses", function(object) {
	# classes.string <- paste(object@classes.type, collapse=", ")
	# cat("Experiment containing ", nrow(object@classes.summary), " samples in ", length(object@classes.type), " different type of classes named: ",classes.string, ". \n \n", sep="")
	# print(object@classes.summary)
# })

# setGeneric("metaData", function(object) standardGeneric("metaData"))
# setMethod("metaData", "MetaboSet", function(object) object@MetaData@Instrumental)
# setGeneric("phenoData", function(object) standardGeneric("phenoData"))
# setMethod("phenoData", "MetaboSet", function(object) object@MetaData@Phenotype)



## Main Software functions:


setDecPar <- function(min.peak.width, min.peak.height=500, noise.threshold=500, avoid.processing.mz=c(73:75,147:149), compression.coef=2, analysis.time=0)
{
	softPar <- new("eRahSoftParameters",algorithm="eRah-OSD", min.peak.width = min.peak.width/60, min.peak.height = min.peak.height, noise.threshold = noise.threshold, avoid.processing.mz = avoid.processing.mz, compression.coef = compression.coef, analysis.time=analysis.time)
	softPar
}

setAlPar <- function(min.spectra.cor, max.time.dist, mz.range=c(70:600))
{
	alPar <- new("eRahAlParameters", algorithm="eRah", min.spectra.cor=min.spectra.cor, max.time.dist=max.time.dist/60, mz.range = mz.range, method="eRah")
	alPar
}

newExp <- function(instrumental, phenotype=NULL, info=character())
{
	#IF es un path:
	#IF path.dir== 
	path.dir <- strsplit(instrumental, split="/")[[1]]
	path.dir <- paste(path.dir[-length(path.dir)], collapse="/")
	
	instrumental.dataframe <- suppressWarnings(try(read.csv(instrumental, sep=";"), silent=T))
	if(class(instrumental.dataframe)=="try-error") stop(attributes(instrumental.dataframe)$condition)
	delete.rows <- apply(instrumental.dataframe,1,function(x) if(x["sampleID"]==""){TRUE}else{FALSE})
	if(any(delete.rows)) instrumental.dataframe <- instrumental.dataframe[-which(delete.rows==T),]
		
	if(is.null(phenotype)) 
	{
		phenotype.dataframe = as.data.frame(NULL)
		warning("No phenotype data have been attached to this experiment.")
	}else{
		phenotype.dataframe <- suppressWarnings(try(read.csv(phenotype, sep=";"), silent=T))
		if(class(phenotype.dataframe)=="try-error") stop(attributes(phenotype.dataframe)$condition)
		## Comprobar almenys que estigui la columna que la relaciona amb la instrumental sampleID
	}
	
	factors.list <- lapply(1:nrow(instrumental.dataframe), function(x){as.data.frame(NULL)})
		
	names(factors.list) <- as.vector(instrumental.dataframe$sampleID)
	ident.list <- as.data.frame(matrix(0,ncol=7, dimnames=list(row=0,col= c("AlignID", "tmean", "Name", "MatchFactor", "CAS", "Formula", "DB.Id"))))
	uni.stats <- as.data.frame(matrix(0,ncol=3, dimnames=list(row=0,col= c("Id", "FoldChangue", "pvalue"))))
	multi.stats <- as.data.frame(matrix(0,ncol=3, dimnames=list(row=0,col= c("Id", "CompoundsInvolved", "pvalue"))))

	al.par <- list()
	id.par <- list()
	soft.par <- list()
	
	stat.parameters <- new("MSResultsParameters", Alignment=al.par, Identification=id.par)
	statistics <- new("Statistics", Univariate = uni.stats, Multivariate = multi.stats)
	MS.Results <- new("Results", Parameters = stat.parameters, Identification = ident.list, Statistics = statistics )
	MS.Data <- new("Data", FeatureList = list(NULL), FactorList = factors.list, Parameters = list(NULL))
	MS.MetaData <- new("MetaData", Instrumental = instrumental.dataframe, Phenotype = phenotype.dataframe, DataDirectory=path.dir)
		
	# Instrumental Slots validation:
	col.correct <- c("sampleID","filename","date","time")
	for(i in 1:length(col.correct))
		if(length(apply(as.matrix(colnames(MS.MetaData@Instrumental)),1,function(x) grep(col.correct[i],x)))==0) stop("Invalid instrumental file. The file must contain at least the following columns: ", paste(col.correct, collapse=", "))

	# Phenotype Slots validation:
	if(!is.null(MS.MetaData@Phenotype))
	{
	col.correct <- c("sampleID","class")
	for(i in 1:length(col.correct))
		if(length(apply(as.matrix(colnames(MS.MetaData@Phenotype)),1,function(x) grep(col.correct[i],x)))==0) stop("Invalid phenotype file. The file must contain at least the following columns: ", paste(col.correct, collapse=", "))
	}
	sample.container <- new("MetaboSet", Info = info, Data = MS.Data, MetaData = MS.MetaData, Results = MS.Results)
	sample.container
}

deconvolveComp <- function(Experiment, decParameters, samples.to.process=NULL)
{
	plotting=FALSE
	Number.of.Samples <- nrow(Experiment@MetaData@Instrumental)
	if(is.null(samples.to.process)) samples.to.process <- 1:Number.of.Samples
	stopifnot(samples.to.process>=1, max(samples.to.process)<=Number.of.Samples, length(samples.to.process)<=Number.of.Samples)
	
	soft.par <- list(min.peak.width = decParameters@min.peak.width, min.peak.height = decParameters@min.peak.height, noise.threshold = decParameters@noise.threshold, avoid.processing.mz = decParameters@avoid.processing.mz,  compression.coef = decParameters@compression.coef, analysis.time = decParameters@analysis.time)
	Experiment@Data@Parameters <- soft.par

	k <- 1
	for(index in samples.to.process)
	{
		cat("\n Deconvolving compounds from",as.character(Experiment@MetaData@Instrumental$filename[index]),"... Processing", k,"/",length(samples.to.process),"\n")  
		Experiment <- processSample(Experiment, index, plotting)
		k <- k + 1
	}
	cat("\n Compounds deconvolved \n")
	Experiment	
}

alignComp <- function(Experiment, alParameters, blocks.size=NULL)
{
	al.par <- list(alignment.algorithm=alParameters@algorithm, min.spectra.cor=alParameters@min.spectra.cor, max.time.dist=alParameters@max.time.dist, mz.range=alParameters@mz.range)
	Experiment@Results@Parameters@Alignment <- al.par
	
	min.spectra.cor <- Experiment@Results@Parameters@Alignment$min.spectra.cor
	max.time.dist <- Experiment@Results@Parameters@Alignment$max.time.dist
	mz.range <- Experiment@Results@Parameters@Alignment$mz.range
	maxMZ <- max(mz.range)
	
	# Experiment@Data@FactorList <- align.factors(Experiment@Data@FactorList, min.spectra.cor, max.time.dist, maxMZ, mz.range)
	# Experiment@Results@Alignment <- create.factorlist.table(Experiment)
	
	if(is.null(blocks.size))
	{
		Experiment@Data@FactorList <- align.factors(Experiment@Data@FactorList, min.spectra.cor, max.time.dist, maxMZ, mz.range)
		Experiment@Results@Alignment <- create.factorlist.table(Experiment)
		#return(Experiment)
	}else{
		
		#blocks.size <- 15
		max.mz <- maxMZ
		Itrt <- length(Experiment@Data@FactorList)/blocks.size
		sequs <- trunc(seq(1, length(Experiment@Data@FactorList), length.out=Itrt))
		sequs[1] <- 0
		
		corresponding.list <- list()
		block.list <- list()
		#i <- 1
	
		for(i in 1:(length(sequs)-1))	
		{
			cat("Aligning block ", i, " of ", length(sequs)-1, "... \n", sep="")
			ghost.object <- Experiment
			ghost.object@Data@FactorList <- Experiment@Data@FactorList[(sequs[i]+1):sequs[(i+1)]]
			factors.list <- ghost.object@Data@FactorList
			ghost.object@Data@FactorList <- align.factors(factors.list, min.spectra.cor, max.time.dist, max.mz, mz.range)
			ghost.factors.list <- create.factorlist.table(ghost.object)
		
			block.list[[i]] <- data.frame(ID=ghost.factors.list$AlignID, RT=ghost.factors.list$tmean, Spectra=ghost.factors.list$Spectra)
			corresponding.list <- c(corresponding.list,lapply(ghost.object@Data@FactorList, function(x) x$AlignID))		
		}
		
		cat("Aligning factors across blocks... \n")
		full.factorlist <- align.factors(block.list, min.spectra.cor, max.time.dist, max.mz, mz.range)
			
		#MaxALID <- max(unlist(lapply(full.factorlist, function(x) x$AlignID)))
		factors.list <- Experiment@Data@FactorList
		if(!(any(unlist(lapply(factors.list,function(x) {is.null(x$AlignID)}))==FALSE)))
		{	
			factors.list <- lapply(factors.list, function(x){
				outp <- cbind(x,matrix(0,nrow=length(x$ID)))
				colnames(outp)[ncol(outp)] <- "AlignID"
				outp
				})
		}else{
			factors.list <- lapply(factors.list, function(x){
				x$AlignID <- rep(0,length(x$ID))
				x
				})
		}	
		
		Experiment@Data@FactorList <- factors.list
			
		free.aligned.slots <- list()		
		for(i in 1:length(full.factorlist))
		{
			for(j in (sequs[i]+1):sequs[(i+1)])
			{
				ID.vct <- sapply(full.factorlist[[i]]$ID, function(x) {x.num <- which(corresponding.list[[j]]==x)
				if(length(x.num)==0) x.num=0
				x.num
				})
				
				#full.factorlist[[i]]$AlignID[which(ID.vct!=0)]	
				#ID.vct[which(ID.vct!=0)]	
			
				Experiment@Data@FactorList[[j]]$AlignID[ID.vct[which(ID.vct!=0)]] <- full.factorlist[[i]]$AlignID[which(ID.vct!=0)]	
				free.aligned.slots[[j]] <- which(full.factorlist[[i]]$AlignID[which(ID.vct!=0)]==0)
			}	
		}
		MaxALID <- max(unlist(lapply(Experiment@Data@FactorList, function(x) x$AlignID)))
		Alid.counter <- MaxALID + 1
		
		for(i in 1:length(free.aligned.slots))
		{
			Experiment@Data@FactorList[[i]]$AlignID[free.aligned.slots[[i]]] <- seq(Alid.counter, Alid.counter + (length(free.aligned.slots[[i]])-1) )
			Alid.counter <- Alid.counter + length(free.aligned.slots[[i]])
		}
		
		cat("Constructing Factor List Table... (This may take a while...)\n")	
		Experiment@Results@Alignment <- create.factorlist.table(Experiment)

	}
	
	Experiment
}


identifyComp <- function(Experiment, id.database=mslib, mz.range=NULL, n.putative=3)
{
	#if(!(any(unlist(lapply(Experiment@Data@FactorList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("Factors must be aligned first")
	
	if(is.null(Experiment@Results@Parameters@Alignment$mz.range) && is.null(mz.range)) stop("A mz.range has to be specified")
	if(is.null(mz.range)) compare.only.mz <- 1:max(Experiment@Results@Parameters@Alignment$mz.range)

	
	id.par <- list(database.name = id.database@name, compare.only.mz = compare.only.mz, n.putative = n.putative)
	Experiment@Results@Parameters@Identification <- id.par
	
	avoid.processing.mz <- Experiment@Data@Parameters$avoid.processing.mz
	maxMZ <- max(compare.only.mz)
	Experiment@Results@Identification <- identify.factors(Experiment, maxMZ, compare.only.mz, avoid.processing.mz, id.database@database, n.putative)
	Experiment
}

processSample <- function(Experiment, index, plotting)
{
	if(Experiment@MetaData@DataDirectory=="") {filename <- as.character(Experiment@MetaData@Instrumental$filename[index])
		}else{filename <- paste(Experiment@MetaData@DataDirectory,"/",Experiment@MetaData@Instrumental$filename[index], sep="")}
	
	sampleObject <- NULL
	sampleObject <- load.file(filename)
	
	# file.extension <- strsplit(as.character(Experiment@MetaData@Instrumental$filename[index]), split="\\.")[[1]]
	# file.type <- file.extension[length(file.extension)]	
	# if(file.type=="cdf") sampleObject <- load.ncdf(filename)
	# if(file.type=="mzXML" || file.type=="xml") sampleObject <- load.xml(filename)
	# if(file.type=="MetaboSet")
	# {
		# load(filename)
		# sampleObject <- new("RawDataParameters", data = sampleRD@data, min.mz = sampleRD@min.mz, max.mz = sampleRD@max.mz, start.time = sampleRD@start.time, mz.resolution = 1)
	# }
		
	Experiment@Data@Parameters$scans.per.second <- sampleObject@scans.per.second
	sampleObject@avoid.processing.mz <- Experiment@Data@Parameters$avoid.processing.mz
	sampleObject@min.peak.width <- Experiment@Data@Parameters$min.peak.width*Experiment@Data@Parameters$scans.per.second*60
	sampleObject@min.peak.height <- Experiment@Data@Parameters$min.peak.height
	sampleObject@noise.threshold <- Experiment@Data@Parameters$noise.threshold
	#sampleObject@moving.window.length <- Experiment@Data@Parameters$moving.window.length*Experiment@Data@Parameters$scans.per.second*60
	#sampleObject@moving.window.overlap <- Experiment@Data@Parameters$moving.window.overlap
	sampleObject@compression.coef <- Experiment@Data@Parameters$compression.coef
	#sampleObject@factor.minimum.sd <- Experiment@Data@Parameters$factor.minimum.sd

	#sampleObject@filter.matrix <- get.filter.matrix(sampleObject)
	
	
		sampleObject <- avoid.processing(sampleObject)
		factor.list <- try(get.factor.list(sampleObject, analysis.window=Experiment@Data@Parameters$analysis.time, plotting), silent=F)
		if(class(factor.list)=="try-error") {factor.list <- as.data.frame(NULL); warning("Unable to extract factors from ", Experiment@MetaData@Instrumental$filename[index], ". Data may be corrupted.", sep="")}
		Experiment@Data@FactorList[[index]] <- factor.list		
	
	
	Experiment
}










