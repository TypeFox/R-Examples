#
# assignLinkageGroups.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Apr, 2013
# first written Oct, 2011
# Contains: assignLinkageGroups, regorganizeMarkersWithin
#

# assignLinkageGroups
#
# DESCRIPTION:
#  Assign linkage groups based on a user supplied known number of chromosomes
#  we can use the genetic map of a cross, or the rf matrix
# OUTPUT:
#  an object of class cross
#
assignLinkageGroups <- function(cross, n.chr, use = c("geno","rf"), ...){
  use 			<- match.arg(use)
  genotype 		<- t(pull.geno(cross))
  genotype[which(is.na(genotype))] <- min(genotype, na.rm = T) - 1 #Why -1?
  
  
  if(use == "geno"){ 
  	clustering 	<- kmeans(genotype, n.chr, nstart=100, ...)
  }else if(use == "rf"){
    cross <- cleanRfs(cross)
    clustering 	<- kmeans(lowerTrng.internal(est.rf(cross)$rf), n.chr, nstart=100, ...)
  }else{ #Added an Else case
  	stop("Please specify if a genetic map or an rf matrix is to be used.")
  }
  reorganizeMarkersWithin(cross, clustering$cluster)
}

lowerTrng.internal <- function(inputMatrix){
  if(ncol(inputMatrix)!=nrow(inputMatrix)) stop("unable to select lower triangle\n")
  outputMatrix <- inputMatrix
  for(i in 1:(nrow(inputMatrix)-1)){
      outputMatrix[i,(i+1):nrow(inputMatrix)] <- outputMatrix[(i+1):nrow(inputMatrix),i]
  }
  invisible(outputMatrix)
}

cleanRfs <- function(cross){
  dataRf 		<- t(est.rf(cross)$rf) #Replaced double for loop with a transform
  rowMinimums	<- apply(dataRf, 1, min, na.rm=T)
  for(i in 1:nrow(dataRf)){
    if(rowMinimums[i] > 0.5){
      cat("dropping", rownames(dataRf)[i], "min:", rowMinimums[i])
      cross <- drop.markers(cross, rownames(dataRf)[i])
    }
  }
  invisible(cross)
}

# reorganizeMarkersWithin
#
# DESCRIPTION:
#  Function to quickly rearrange all the markers in a cross based on any ordering vector
# OUTPUT:
#  an object of class cross
#
reorganizeMarkersWithin <- function(cross, ordering){
  cross 			<- clean(cross) #What is this clean function? Where is it?
  chrType 			<- rep(sapply(cross$geno, class), sum(nmar(cross)))
  crossType 		<- class(cross)[1]
  genotype 			<- pull.geno(cross)
  newChromosomes 	<- sort(unique(ordering))
  n.newChromosomes 	<- length(newChromosomes)
  
  cross$geno 		<- vector("list", n.newChromosomes)
  names(cross$geno) <- newChromosomes
  
  for (i in 1:n.newChromosomes) {
    selectedMarkers <- which(ordering == newChromosomes[i])
    cross$geno[[i]]$data <- genotype[, selectedMarkers, drop = FALSE]
    cross$geno[[i]]$map <- seq(0, by = 10, length = length(selectedMarkers))
      
    if (crossType == "4way") stop("4way is not supported by pheno2geno.")
      	
    names(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
    uniqueChrType <- unique(chrType[which(ordering == newChromosomes[i])])
      
    if (length(uniqueChrType) > 1){          
      warning("Problem with linkage group ", i, ": A or X?\n", paste(uniqueChrType, collapse = " "))
    }else{ 
      class(cross$geno[[i]]) <- uniqueChrType
    }
  }	
 
  cross <- est.rf(cross)
  return(cross)
}
