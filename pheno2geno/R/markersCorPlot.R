#
# markersCorPlot.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Mar, 2011
# Contains: markersCorPlot, ascendingMaptoJigSawMap, getChrOffsets.internal
#           getMarkerOffsets, getMarkerOffsetsFromMap, chromCorMatrix
#           getPopulationOffsets.internal 
#

# markersCorPlot
#
# DESCRIPTION:
#  function to create new map and save it in cross object
# PARAMETERS:
#   - population - object of class population
#   - orde - object of class population
#   - n.chr - expected number of linkage groups
#   - use - expected number of linkage groups
#   - verbose - be verbose
# OUTPUT:
#  An object of class cross
#
markersCorPlot <- function(cross, population, map=c("genetic","physical"), cmBetween=25, comparisonMethod = c(sumMajorityCorrelation,majorityCorrelation,meanCorrelation,majorityOfMarkers), chr, show.legend=FALSE, verbose=TRUE){
  if(missing(cross)) print("No cross object provided.")
  if(missing(population)) print("No population object provided.")
  map <- match.arg(map)  # Check
  if(map=="genetic"){
    originalMap <- population$maps$genetic
  }else{
    originalMap <- population$maps$physical
  }
  comparisonMethod <- defaultCheck.internal(comparisonMethod,"comparisonMethod",4,sumMajorityCorrelation)
  if(is.null(originalMap)) stop("no ",map," map provided!")  
  
  ### getting offsets for each chromosome on both maps
  offsets1 <- getPopulationOffsets.internal(population,originalMap,cmBetween)
  n.originalChrom <- length(offsets1)-1
  offsets2 <- getChrOffsets.internal(cross,cmBetween)
  n.newChrom <- length(offsets2)-1
  
  if(n.originalChrom<n.newChrom){
	offsets1 <- c(offsets1,rep(0,(n.newChrom- n.originalChrom)))
  }else if( n.originalChrom>n.newChrom){
	offsets2 <- c(offsets2,rep(0,( n.originalChrom-n.newChrom)))
  }
    
  ### global offsets
  global_offset <- NULL
  for(x in 1:length(offsets1)){
    global_offset <- c(global_offset,max(offsets1[x],offsets2[x]))
  }
  
  ### positions of markers (absolute - with offsets)
  mloc_original <- getMarkerOffsets(cross,global_offset[1:n.newChrom],cmBetween)
  mloc_o <- getMarkerOffsetsFromMap(originalMap,global_offset[1:n.originalChrom],cmBetween)

  ### limits of plot

  
  ### summary offsets
  sum_gl_off <- NULL
  for(x in 1:length(global_offset)){
    sum_gl_off <- c(sum_gl_off,sum(global_offset[1:x]))
  }
  if(missing(chr)){
    m_max <- max(mloc_o,mloc_original)
    m_min <- min(mloc_o,mloc_original)
  }else{
    m_max <- sum_gl_off[max(chr)+1]
    m_min <- sum_gl_off[min(chr)]
  }
  ### preparing chrom to chrom cor matrix for use in the background
  genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross, population, FALSE)
  chromToChromMatrix <- comparisonMethod(cross,originalMap,population)
  maximum <- max(chromToChromMatrix)
  
  ### setting plot canvas
  plot(c(m_min,m_max),c(m_min,m_max),type='n',xlab="Original map",ylab="New map",main="Comparison of genetic maps", xaxt="n", yaxt="n")
  ### background
  for(i in 1:(n.originalChrom)){
    for(j in 1:(n.newChrom)){
        cur_col <- (maximum-(chromToChromMatrix[i,j]))/maximum
        rect(sum_gl_off[i],sum_gl_off[j],sum_gl_off[i+1],sum_gl_off[j+1],lty=0,col=rgb(cur_col,cur_col,cur_col))
    }
  }
  size_ <- (m_max-m_min)
  ### markers on new map
  points(cbind(mloc_o,mloc_o+0.01*size_),pch=20,col="green",cex=1.5,lwd=2)
  ### markers on original map
  points(cbind(mloc_original+0.01*size_,mloc_original),pch=20,col="red",cex=1.5,lwd=2)
  ### gris
  abline(v=sum_gl_off[-length(sum_gl_off)],lty=2)
  abline(h=sum_gl_off[-length(sum_gl_off)],lty=2)
  ### chromosome labels and tics
  labelsPos <- vector(mode="numeric",length(sum_gl_off)-1)
	for(i in 1:length(sum_gl_off)-1){
		labelsPos[i] <- (sum_gl_off[i] + sum_gl_off[i+1])/2
	}
	axis(1, at = sum_gl_off[-length(sum_gl_off)],labels = FALSE)
	axis(1, at = labelsPos[1:n.originalChrom],labels = names(table(originalMap[,1])), lwd = 0, tick = FALSE)
	axis(2, at = sum_gl_off[-length(sum_gl_off)],labels = FALSE)
	axis(2, at = labelsPos[1:n.newChrom],labels = chrnames(cross), lwd = 0, tick = FALSE)
  if(show.legend){
    maximum <- round(maximum) 
    text_ <- c(0,maximum*0.2,maximum*0.4,maximum*0.6,maximum*0.8,maximum)
    colors_ <- c(rgb(1,1,1),rgb(0.8,0.8,0.8),rgb(0.6,0.6,0.6),rgb(0.4,0.4,0.4),rgb(0.2,0.2,0.2),rgb(0,0,0))
    legend(c(m_min,size_*0.3),c(size_*0.67,size_),text_,fill=colors_,bg="white",cex=0.8)
    legend(c(m_min,size_*0.3),c(size_*0.67,size_*0.77),c("original","new"),col=c("red","green"),cex=0.8,pch=20,bty="n",lwd=4)
  }
  invisible(chromToChromMatrix)
}

############################################################################################################
#									*** markersCorPlot ***
#
# DESCRIPTION:
# 	function to create new map and save it in cross object
# 
# PARAMETERS:
# 	population - object of class population
# 	orde - object of class population
# 	n.chr - expected number of linkage groups
# 	use - expected number of linkage groups
#	verbose - be verbose
#
# OUTPUT:
#	an object of class cross
#
#
############################################################################################################
getChrOffsets.internal <- function(cross, cmBetween){
  offsets <- unlist(lapply(pull.map(cross),max))
  offsets <- offsets+cmBetween
  offsets <-c(0,offsets)
  offsets
}


#From IQTL by Danny Arends, SHOULD NOT MODIFY
getMarkerOffsets <- function(cross, offsets, cmBetween=25){
  if(missing(offsets))offsets <- getChrOffsets.internal(cross,cmBetween)
  cnt <- 1
  myoffsets <- NULL
  for(x in nmar(cross)){
    myoffsets <- c(myoffsets,rep(sum(offsets[1:cnt]),x))
    cnt <- cnt+1
  }

  mlocations <- myoffsets + as.numeric(unlist(pull.map(cross)))
  mlocations
}

############################################################################################################
#									*** markersCorPlot ***
#
# DESCRIPTION:
# 	function to create new map and save it in cross object
# 
# PARAMETERS:
# 	population - object of class population
# 	orde - object of class population
# 	n.chr - expected number of linkage groups
# 	use - expected number of linkage groups
#	verbose - be verbose
#
# OUTPUT:
#	an object of class cross
#
#
############################################################################################################
getMarkerOffsetsFromMap <- function(map, offsets, cmBetween=25){
  cnt <- 1
  myoffsets <- NULL
  for(x in table(map[,1])){
    myoffsets <- c(myoffsets,rep(sum(offsets[1:cnt]),x))
    cnt <- cnt+1
  }

  mlocations <- myoffsets + as.numeric(map[,2])
  mlocations
}


############################################################################################################
#									*** markersCorPlot ***
#
# DESCRIPTION:
# 	function to create new map and save it in cross object
# 
# PARAMETERS:
# 	population - object of class population
# 	orde - object of class population
# 	n.chr - expected number of linkage groups
# 	use - expected number of linkage groups
#	verbose - be verbose
#
# OUTPUT:
#	an object of class cross
#
#
############################################################################################################
ascendingMaptoJigSawMap <- function(mapToProcess,verbose=FALSE){
  if(is.null(mapToProcess)) stop("No map provided!")
  for(x in unique(mapToProcess[,1])){
    if(verbose) cat("Processing chromosome:",x,"\n")
    offsetOfCurrentChromosome <- min(mapToProcess[which(mapToProcess[,1]==x),2])
    mapToProcess[which(mapToProcess[,1]==x),2] <- mapToProcess[which(mapToProcess[,1]==x),2]-offsetOfCurrentChromosome
  }
  invisible(mapToProcess)
}

############################################################################################################
#									*** markersCorPlot ***
#
# DESCRIPTION:
# 	function to create new map and save it in cross object
# 
# PARAMETERS:
# 	population - object of class population
# 	orde - object of class population
# 	n.chr - expected number of linkage groups
# 	use - expected number of linkage groups
#	verbose - be verbose
#
# OUTPUT:
#	an object of class cross
#
#
############################################################################################################
getPopulationOffsets.internal <- function(population, originalMap, cmBetween){
  minima <- NULL
  for(x in unique(originalMap[,1])){
    minima <- c(minima,max(originalMap[which(originalMap[,1]==x),2]))
  }
  minima <- minima + cmBetween
  minima <- c(0,minima)
  invisible(minima)
}
