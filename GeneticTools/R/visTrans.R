visTrans <- function(snpGene, geneAnnot=NULL){
 
 if(!is.null(geneAnnot)) annotTrack <- TRUE

# Extract the genes
  genes <- unique(as.character(snpGene[,5]))

# Get the Annotations from it 
  annot <- geneAnnot[is.element(geneAnnot[,1],genes),]

# Input checks and cleaning
  if(length(genes)!=nrow(annot)) warning("Not all genes have been annotated!\n")
  removeThese <- which((is.na(snpGene[,2]))==TRUE)
  if(length(removeThese)>0) snpGene <- snpGene[-removeThese,]
  
# Determine constants for the plot
  temp <- c(snpGene[,3],annot[,3],annot[,4])
  xlim <- c(min(temp),max(temp))  / 10^6
  rcolors <- rainbow(1000,end=4/6)
 
# Determine the color of each gene (depending on the position on the chromosome)
  geneColors <- rcolors[floor((apply(cbind(annot[,3],annot[,4]),1,mean)/10^6 - xlim[1])/ (xlim[2]-xlim[1]) * length(rcolors))]
 
  snpGene[,3] <- snpGene[,3] / 10^6
  annot[,3] <- annot[,3] / 10^6
  annot[,4] <- annot[,4] / 10^6
  
# Create the layout for the plot
  if(annotTrack==TRUE)
  {
    layout(matrix(c(1,2),ncol=1),heights=c(1,2))
    par(mar=c(0,4,0,2))
  }
  
# Create the upper plot
  if(annotTrack==TRUE)
  {
    plot(-100,-100,xlim=xlim, ylim=c(0,1), xlab="", xaxt="n", yaxt="n",ylab="")
  } else {
    plot(-100,-100,xlim=xlim, ylim=c(0,1), xlab=paste("Chromosomal Position"),yaxt="n",ylab="")
  }

# Now connect the genes and the SNP 
  for(i in 1:nrow(snpGene)){
    temp <- snpGene[i,]
  # Get the temporabry locations
    thisGene <- annot[,1]==as.character(temp[,5])
    geneLoc <- mean(as.numeric(annot[thisGene,c(3,4)]))
    snpLoc <- snpGene[i,3]
    ellColor <- geneColors[thisGene]
  # Draw the ellipse 
    gsDist <- abs(geneLoc-snpLoc)
    draw.ellipse(mean(c(geneLoc,snpLoc)),0,gsDist/2,gsDist/(xlim[2]-xlim[1]),segment=c(0,180),border=ellColor)
  }

# Add the genes 
#   for(i in 1:length(genes))
#   {
#     rect(annot[i,3],0,annot[i,4],0.1,col="black")
#   }
# Add the SNPs 
#   for(i in 1:nrow(snpGene))
#   {
#     rect(snpGene[i,3],0,snpGene[i,3],0.05,col="red",border="red")
#   }

  colnames(annot) <- c("Name","Chr","Start","Stop")
  annot[,3] <- annot[,3] * 10^6
  annot[,4] <- annot[,4] * 10^6
  annot <- annot[,c(2,3,4,1)]
# Now add the lower plot  
 if(annotTrack==TRUE)
 {
     par(mar=c(5,4,0,2))
     plotGeneLocation(annot,chrInfo=2,xlim=xlim,minDistance=1.5)
 }
}