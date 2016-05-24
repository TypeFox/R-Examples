#
# saveGff.r
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Nov, 2012
# first written Nov, 2012
# Contains: saveGff, findrecombinations
#

#  saveGff
#
# DESCRIPTION:
#  Saves a gff file containing physical locations for markers for use in genome viewers.
# PARAMETERS:
#   - population -  An object of class cross.
#   - population -  An object of class population.
#   - gffFileCore - an output file
#   - verbose - Be verbose
# OUTPUT:
#  NONE
#
save.gff <- function(cross, map.physical, ind, gffFileCore="population", verbose=FALSE){
  if(missing(ind)) ind <- 1:nind(cross)
  if(missing(map.physical)) stop("Provide map.physical.")
  if(missing(cross)) stop("Provide an object of class cross.")
  markersCross <- markernames(cross)
  markers <- markersCross[which(markersCross %in% rownames(map.physical))]
  if(is.integer0(markers) || is.null(markers)) stop("No physical locations for any of the markers in the cross object!\n")
  #crossNr <- drop.dupmarkers(cross,verbose=FALSE)
  #markersCrossNr <- markernames(crossNr)
  #markersNr <- markersCrossNr[which(markersCrossNr %in% rownames(population$maps$physical))]
  #if(is.integer0(markersNr) || is.null(markersNr)) stop("No physical locations for any of the non redundant markers in the cross object!\n")
  genos <- pull.geno(cross)
  for(x in ind){
    filename <- paste(gffFileCore,"_ind",x,".gff",sep="")
    if(verbose) cat("===",filename,"===\n")
    cat("##gff-version 3\n",file=filename)
    haplotypes <- makeHaplos.internal(map.physical,genos[x,])
    for(marker in 1:nrow(haplotypes)){
          cat("Chr",as.numeric(haplotypes[marker,2]),"\t.\tmarker-",haplotypes[marker,1],"\t",as.numeric(haplotypes[marker,3]),"\t",as.numeric(haplotypes[marker,4]),"\t100\t+\t.\tID=",marker,"\n",file=filename,append=TRUE,sep='')
    }
  #if(verbose) cat("Saved ",length(markersNr),"non-redundant and",length(markers)-length(markersNr),"redundant markers. In total:",length(markers),"markers.\n")
}
}

makeHaplos.internal <- function(map,genotype){
  blues <- matrix(c(1,0,0),1,3)
  reds <- matrix(c(1,0,0),1,3)
  blueRow <- 1
  redRow <- 1
  map <- map[order(map[,1]),]
  firstVal <- genotype[1]
  for(chr in unique(map[,1])){
    chrMap <- map[which(map[,1]==chr),]
    chrMap <- chrMap[order(chrMap[,2]),]
    markers <- rownames(chrMap)
    markerNr <- 1
    oldMarkerVal <- genotype[markers[1]]
    while(markerNr < length(markers)){
      marker <- markers[markerNr]
      markerVal <- genotype[marker]
      if(markerVal!=oldMarkerVal){
        if(markerVal==firstVal){
          reds[redRow,3] <- chrMap[marker,2]
          blues <- rbind(blues,c(chr,chrMap[marker,2],chrMap[marker,2]))
          blueRow <- blueRow+1
        }else{
          blues[blueRow,3] <- chrMap[marker,2]
          reds <- rbind(reds,c(chr,chrMap[marker,2],chrMap[marker,2]))
          redRow <- redRow+1
        }
      }
      oldMarkerVal <- genotype[marker]
      markerNr <- markerNr+1
    }
    if(markerVal==firstVal){
      blues[blueRow,3] <- chrMap[marker,2]
      if(chr<max(unique(map[,1]))){
        chrMap <- map[which(map[,1]==chr+1),]
        chrMap <- chrMap[order(chrMap[,2]),]
        chrNewVal <- genotype[rownames(chrMap)[1]]
        if(chrNewVal==markerVal){
          blues <- rbind(blues,c(chr+1,0,0))
          blueRow <- blueRow+1
        }else{
          reds <- rbind(reds,c(chr+1,0,0))
          redRow <- redRow+1
        }
      }
    }else{
      reds[redRow,3] <- chrMap[marker,2]
      if(chr<max(unique(map[,1]))){
        chrMap <- map[which(map[,1]==chr+1),]
        chrMap <- chrMap[order(chrMap[,2]),]
        chrNewVal <- genotype[rownames(chrMap)[1]]
        chrNewVal <- genotype[rownames(map[which(map[,1]==chr+1),])[1]]
        if(chrNewVal==markerVal){
          reds <- rbind(reds,c(chr+1,0,0))
          redRow <- redRow+1
        }else{
          blues <- rbind(blues,c(chr+1,0,0))
          blueRow <- blueRow+1
        }
      }
    }
  }
  allres <- rbind(cbind(rep("blue",nrow(blues)),blues),cbind(rep("red",nrow(reds[-1,])),reds[-1,]))
  allres <- allres[order(as.numeric(allres[,2])),]
  output <- NULL
  for(chr in unique(map[,1])){
    cur <- allres[which(allres[,2]==chr),]
    cur <- cur[order(as.numeric(cur[,3])),]
    output <- rbind(output,cur)
  }
  invisible(output)
}