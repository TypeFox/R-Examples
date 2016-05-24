plotGeneLocation <- function(bed,chrInfo=NULL,xlim=NULL,order=NULL,mc.cores=6,annot=TRUE, minDistance=NULL){ 
   
  # In case the user gives a xlim object use this information, otherwise plot the whole chromosome or if that is not available,
  # plot the region given from the bed file
  # No xlim object given:
  if(is.null(xlim)){
    # No chr info given
    if(is.null(chrInfo)){
      xlim <- c(min(as.numeric(bed[,2]))/10^6-0.5,max(as.numeric(bed[,3]))/10^6+0.5)
    } else {
      xlim <- c(0,chrInfo/10^6)
    }
  }
  # Extract all starting and stoping positions in MB
  start <- as.numeric(bed[,2])/10^6
  stop <- as.numeric(bed[,3])/10^6
  
  # Reduce the bed file to the area of interest:
  bed <- bed[(start>=xlim[1]) & (stop<=xlim[2]),]
  
  bed <- bed[order(bed[,4]),]
  
  # Update to correct values
  start <- as.numeric(bed[,2])/10^6
  stop <- as.numeric(bed[,3])/10^6
  
  # Now visualize all genes (by horizontal lines)
  createGeneBoundaries <- function(bed,gene){
    exons <- bed[bed[,4]==gene,c(2,3)]
    
    # Bring it in the right format:
    exons <- as.numeric(as.vector(as.matrix(exons)))
    return(data.frame(Start=min(exons), Stop=max(exons), Gene=as.character(gene)))
  }
  
  # Get the Gene locations:
  tempGenes <- as.character(unique(bed[,4]))
  if(is.null(order)){
    geneLocations <- mclapply(tempGenes, createGeneBoundaries, bed=bed,mc.cores=mc.cores)
    bind.ith.rows <- function(i) do.call(rbind, lapply(geneLocations, "[", i, TRUE))
    nr <- nrow(geneLocations[[1]])
    geneLocationsJoined <- lapply(1:nr, bind.ith.rows)[[1]]
  } else {
    geneLocationsJoined <- bed[,2:4]
  }
  
  StartGenes <- as.numeric(geneLocationsJoined[,1])/10^6
  StopGenes <- as.numeric(geneLocationsJoined[,2])/10^6

  # Now we have to calculate the right layer for each gene, 
  # The minimum distance for two genes on the same layer should be 10% of the total distance:
  if(is.null(minDistance)) minDistance <- (xlim[2]-xlim[1])/20
  
  if(is.null(order)){
    # Vector of layer positions for each gene:
    layer <- rep(1,length(StartGenes))
    # Subfunction that increase the layer for the current starting gene
    checkLayer <- function(pos,start){
      for(i in start:nrow(pos)){
	# If the minimum distance of two following values isn't at least the min distance, increase the layer for the next 
	if(((pos[i,1]-pos[start-1,2]) < minDistance) && (layer[i]==layer[start-1])) layer[i] <<- layer[i] + 1
      }
    }
    temp <- cbind(StartGenes,StopGenes)
    for(i in 2:(length(StartGenes))) checkLayer(temp,i)
    
    getMult <- function(x){
      possibilities <- length(unique(x))
      mult <- rep(1,possibilities)
      mPos <- 1
      for(i in 2:length(x)) ifelse((x[i]==x[i-1]), mult[mPos] <- mult[mPos]+1 , mPos <- mPos + 1)
      mult
    }
    
    multiplicities <- getMult(as.character(bed[,4]))
    exonLayer <- rep(layer,multiplicities)
  } else {
    exonLayer <- rep(0,nrow(geneLocationsJoined))
    for(i in 1:length(order)){
      layer <- 1:length(order)
      layerPos <- which((is.element(as.character(geneLocationsJoined[,3]),order[i])==TRUE))
      exonLayer[layerPos] <- layer[i]
    }
  }
  # create the basic plot 
  plot(-100,-100,xlim=xlim,ylim=c(0,max(layer)+1),yaxt="n",xlab="Chromosomal Position",ylab="")
  
  # Plot all exons and genes
  if(is.null(order))rect(StartGenes,layer-0.025,StopGenes,layer+0.025,col="grey",border="grey")
  rect(start,exonLayer-0.1,stop, exonLayer+0.1,col="black")
  
  # Add the Gene names
  if(annot==TRUE){
     text(StartGenes,layer+0.3,unique(as.character(bed[,4])))
  } else {
    if(!is.null(order)){
      axis(2,at=1:length(order),labels=order,las=2)
      for(i in 1:(length(order)-1)){
	lines(c(xlim[1]-10,xlim[2]+10),c(i+0.5,i+0.5),lty="dotted")
      }
    }
  }
  
  # Add still the labels of the genes and then also the different layers (this might be the tricky part...)
}

