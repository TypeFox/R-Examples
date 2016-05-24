
#  ParentOffspring.R

ParentOffspring = function(genoMat, plot=TRUE){

  #############################################################################
  #
  #  Conduct the parent-offspring test using monomorphic SNP markers
  #
  #  Input: 
  #
  #       genoMat: matrix or data frame with n rows and m columns where the rows are for the SNP markers, 
  #           and the first two columns are for the parents and the remaining columns are for the offspring. 
  #           Each genoMat value is a string of length 2. For example,
  #           "GG" "GG" "GG" "CG" "CG" "GG"
  #           "AA" "AA" "AA" "--" "AA" "AA"
  #           "AA" "AA" "GG" "AA" "AA" "AA"
  #           "GG" "AA" "GG" "GG" "GG" "GG"  
  #           Note the procedure will use only the rows with the monomorphic parents, i.e., 
  #           some rows such as the last one here will be ignored.
  #
  #       plot:  logical value for ploting the similarity graph. The default value is TRUE. 
  #
  #  Return: 
  #       similarity: vector of similarity for each offspring
  #
  #########################################################################
  
  # convert a data frame into a matrix, with all values as strings instead of factors
  if(class(genoMat)=="data.frame")  genoMat=as.matrix(genoMat)
  
  #check input errors
  dimension = dim(genoMat)
  if(dimension[1]<1) stop("the input genoMat is empty")
  if(dimension[2]<3) stop("the input genoMat should have at least 3 columns")
  aa=unlist(lapply(unlist(genoMat), is.character))
  if(min(aa)==0) stop("the input genoMat has noncharacter values")
  if(min(nchar(unlist(genoMat))==2) ==0) stop("the input genoMat has values which are not 2-letter strings")  
    
  # Extract all nonmissing monomorphic markers
  geno <- genoMat
  # colnames(geno)=c("P1", "P2")
  # Keep only the rows where P1 (parent 1) = P2 (parent 2)
  geno <- geno[(geno[,1]) == (geno[,2]), ]
  # Keep only the rows where P1 has two identical letters (biallelic fashion)
  geno <- geno[(substr(geno[,1],1,1)) == (substr(geno[,1],2,2)),]
  # Remove the rows were P1 are missing
  geno <- geno[(geno[,1])!="--",]
  dimension <- dim(geno)
  markersN <- dimension[1]
  offspringN <- dimension[2] - 2
  
  #check if  markersN is positive
  if(markersN < 1) stop("No monomorphic markers can be used")
  
  # Compute the similarity matrix
  parent <- geno[,1]
  offspring <- geno[, -(1:2)]
  markersValidN <- rep(0, offspringN) # for each offspring
  similarityMat <- matrix(0, markersN, offspringN)
  for (i in 1:offspringN){
    # Calculate the numerical (0,1,2) vector for each offspring
    similarityMat[, i] <- (substr(parent,1,1) == substr(offspring[, i], 1, 1)) +
      (substr(parent, 1, 1) == substr(offspring[, i], 2, 2))
    
    # Check the missing genoMat values
    for (j in 1:markersN){
      if ((substr(offspring[j, i], 1, 1) =="-") ||(substr(offspring[j, i], 2, 2) == "-") ){
        similarityMat[j, i] <- NA
      }  
    }
    
    # Calculate the number of valid markers for this offspring
    markersValidN[i] <-  markersN - sum(is.na(similarityMat[, i]))
  }
  similarityMat <- similarityMat/2;
  similarity <- round(colMeans(similarityMat, na.rm=T),2)
  
  # Plot the sorted similarity 
  # pdf("sorted.similarityall-13444.pdf", paper="USr")
  if(plot==TRUE){
    plot(1:offspringN, sort(similarity, decreasing=T), pch='*',cex=2,
         col='red', ylim=c(0.5,1), xlab="offspring sorted by similarity", 
         ylab='similarity', axes=FALSE)
    axis(side = 1, at = 1:offspringN)
    axis(side = 2, at = c(round((5:10)/10, 1), 0.95, 0.99))
    #Add dd dashed blue horizontal lines at y=0.9, 0.95 and 0.99
    abline(h=c(0.9, 0.95,0.99),lty=2,col="blue", lwd=2)
    title("Similarity to the Parents")
    box()
    #dev.off()
  }
  return(similarity)
  
}