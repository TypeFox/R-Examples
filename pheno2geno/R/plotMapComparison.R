#
# plotMapComparison.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Nov, 2011
# Contains: plotMapComparison, getChromosome.internal, getYLocs.internal
#           makeChromPal.internal, chromosomesLengths.internal
#

# plotMapComparison
#
# DESCRIPTION:
#  Boxplot of data for selected markers + points of founders mean for each marker
# PARAMETERS:
#   - cross - object of R/qtl cross type
#   - coloringMode - 1 - rainbow colors 2 - black for cis and red for trans located markers
#   - map - which map should be used for comparison:
#     - genetic - genetic map from cross$maps$genetic
#     - physical - physical map from cross$maps$physical
# OUTPUT:
#  Plot
#
plotMapComparison <- function(cross, population, map=c("genetic","physical"), chr){
  #*******objects containing all information needen for function execution*******
  map <- defaultCheck.internal(map,"map",2,"genetic")
	if(map=="genetic"){
		cur_map <- population$maps$genetic
	}else{
		cur_map <- population$maps$physical
	}
	if(is.null(cur_map)) stop("no ",map," map provided!")
  if(!missing(chr)){
    if(any(!(chr%in%table(cur_map[,1])))) stop ("Wrong subset of chromosomes selected, doesn't match chromosomes from original map.")
    if(any(!(chr%in%1:nchr(cross)))) stop ("Wrong subset of chromosomes selected, doesn't match chromosomes from new map.")
  }
	ys <- getYLocs.internal(cross,chr)
	ys[[1]] <- mapMarkers.internal(ys[[1]],cur_map,1)
	xs <- mapMarkers.internal(cur_map,ys[[1]],1)
	#*******chromosomes lengths*******
	referenceChrom <- chromosomesLengths.internal(xs)
	xs[,2] <- xs[,2] + referenceChrom[xs[,1]]
	
	#*******positions of markers*******
	predictedLocs <- ys[[1]][,-1]
	referenceLocs <- xs[,2]
	
	#*******chromosomes lengths*******
	predictedChrom <- ys[[2]]
	
	#*******chromosome labels*******
	predictedChromLabels <- names(table(ys[[1]][,1]))
	referenceChromLabels <- names(table(xs[,1]))
	
	#*******chromosome labels positions*******
	predictedChromPos <- vector(mode="numeric",length(predictedChrom)-1)
	for(i in 1:length(predictedChrom)-1){
		predictedChromPos[i] <- (predictedChrom[i] + predictedChrom[i+1])/2
	}
	predictedChromPos[length(predictedChrom)] <- (predictedChrom[length(predictedChrom)] + max(predictedLocs))/2
	
	referenceChromPos <- vector(mode="numeric",length(referenceChrom)-1)
	for(i in 1:length(referenceChrom)-1){
		referenceChromPos[i] <- (referenceChrom[i] + referenceChrom[i+1])/2
	}
	#referenceChromPos[length(referenceChrom)] <- (referenceChrom[length(referenceChrom)] + max(referenceLocs))/2
	#*******color palette*******
    color <- makeChromPal.internal(ys[[1]],xs)
	#*******plotting points*******
	plot(x=referenceLocs, y=predictedLocs, xlim=c(min(referenceLocs),max(referenceLocs)), ylim=c(min(predictedLocs),max(predictedLocs)),
		xaxt="n", yaxt="n", col=color[[1]], pch=color[[2]], cex=1.5, xlab="Reference map", ylab="Predicted map", main="Comparison of genetic maps")
	#*******adding chromosome labels and tics*******
	if(!missing(chr)){
		referenceChrom <- referenceChrom[chr]
		referenceChromPos <- referenceChromPos[chr]
		predictedChrom <- predictedChrom[chr]
		predictedChromPos <- predictedChromPos[chr]
	}
	axis(1, at = referenceChrom[-1],labels = FALSE)
	axis(1, at = referenceChromPos,labels = referenceChromLabels, lwd = 0, tick = FALSE)
	axis(2, at = predictedChrom[-1],labels = FALSE)
	axis(2, at = predictedChromPos,labels = predictedChromLabels, lwd = 0, tick = FALSE)
	
	#*******adding marker tics*******
	axis(1, at = referenceLocs,labels = FALSE)
	axis(2, at = predictedLocs,labels = FALSE)
	#*******adding lines marking chromosomes ends*******
	for(x in 2:length(referenceChrom)){
		abline(v=sum(referenceChrom[x]),lty=2)
	}
	for(y in 2:length(predictedChrom)){
		abline(h=predictedChrom[y],lty=2)
	}
}

############################################################################################################
#									*** getChromosome.internal ***
#
# DESCRIPTION:
#	subfunction of plotMapComparison, returning list of chromosome numbers for all markers in cross
# 
# PARAMETERS:
# 	cross - object of class cross 
# 
# OUTPUT:
#	vector of numbers
#
############################################################################################################
getChromosome.internal <- function(cross,chr){
	if(missing(chr)){
		invisible(rep(1:length(nmar(cross)),nmar(cross)))
	}else{
		invisible(rep(chr,nmar(cross)[chr]))
	}
}

############################################################################################################
#									*** getYLocs.internal***
#
# DESCRIPTION:
#	subfunction of plotMapComparison, returning list of location of all markers in cross
# 
# PARAMETERS:
# 	cross - object of class cross 
#
# OUTPUT:
#	vector of numbers
#
############################################################################################################
getYLocs.internal <- function(cross,chr){
	locs <- lapply(cross$geno,function(x){as.numeric(x$map)})
	chrlength <- as.vector(unlist(lapply(locs,max)),mode="numeric")
	locs <- as.numeric(unlist(locs))
	summaryLengths <- rep(0,length(chrlength))
	if(length(chrlength)>1){
		for(x in 2:length(chrlength)){
			summaryLengths[x] <- max(chrlength[x-1]) + summaryLengths[x-1] + 0.15 * max((chrlength))
		}
	}
	chrids <- getChromosome.internal(cross)
	result <- matrix(0,length(locs),2)
	result[,1] <- chrids
	result[,2] <- summaryLengths[chrids]+locs
	rownames(result) <- markernames(cross)
	if(!missing(chr)) result <- result[which(result[,1]%in%chr),]
	invisible(list(result,summaryLengths))
}

############################################################################################################
#									*** makeChromPal.internal ***
#
# DESCRIPTION:
#	subfunction of plotMapComparison, returning color pallete (rainbow colors)
# 
# PARAMETERS:
# 	ys1 - object of plotMapComparison, containing info about predicted map
# 	xs - object of plotMapComparison, containing info about reference map
#
# OUTPUT:
#	list of vector of colors (characters) and vector of numbers (symbol identifiers)
#
############################################################################################################
makeChromPal.internal <- function(ys1,xs){
	color <- vector(mode="character",nrow(ys1))
	names(color) <- rownames(ys1)
	symbol <- vector(mode="numeric",nrow(xs))
	names(symbol) <- rownames(xs)
	#cl <- topo.colors(length(table(ys1[,1])))
  cl <- c("red","green","blue")
	for(i in rownames(ys1)){
		color[i] <- cl[ys1[i,1]%%3+1]
		symbol[i] <- 19
	}
	invisible(list(color, symbol))
}

############################################################################################################
#									*** chromosomesLengths.internal ***
#
# DESCRIPTION:
#	function calculating lengths of chromosomes in given map
# 
# PARAMETERS:
# 	map - genetic or physical map (matrix with two cols - 1 - chromome nr, 2- postion on chromosome),
#		rownames - names of markers
#
# OUTPUT:
#	vector of lengths of chromosomes
#
############################################################################################################
chromosomesLengths.internal <- function(map){
	lengths <- vector(mode="numeric",length=(max(unique(map[,1]))+1))
	lengths[1] <- 0
	for(i in unique(map[,1])){
		lengths[i+1] <- max(map[which(map[,1]==i),2]) + lengths[i]
	}
	invisible(lengths)
}
