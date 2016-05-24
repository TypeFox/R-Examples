updateGeneTable <-
function(genetable){
  chromosomerows <- grep("CHR", genetable[,1], ignore.case=T, value=F)
  # The rows (indices) at where source (chromosome) information is located
  chromosomeheaders <- genetable[chromosomerows,]
  # The actual rows
  totalrows <- nrow(genetable)
  resultslist <- list()
  
  for(chr in 1:length(chromosomerows)){ # For each chromosome
    if(chr == length(chromosomerows)){ # If this is the last chromosome in the file
      goto <- totalrows # Include genes down to the last one
    }
    else{
      goto <- chromosomerows[(chr+1)]-1 # Include genes down to the next chromosome 
    }
    chrlength <- genetable[chromosomerows[chr],5]
    this_genes <- genetable[(chromosomerows[chr]+1):goto,]

    new_this_genes <- this_genes
    # Insert the first intergenic region if chromosome does not have a gene at position 1
    k <- 1
    if(!this_genes$Left[1]==1){
      new_this_genes <- insertRow(new_this_genes,newrow=c("IG","",1,1,this_genes$Left[1]-1,this_genes$Left[1]-1),index=1)
      k <- k + 1
    }
    # Insert the intergenic regions between genes throughout the chromosome
    for(i in 1:nrow(this_genes)){
      if(i == nrow(this_genes)) break # If we are at the last line
      if(this_genes$Left[i+1] - this_genes$Right[i] > 1){ # If there is an intergenic space between two genes
        k <- k+1
        new_this_genes <- insertRow(new_this_genes,newrow=c("IG","",1,this_genes$Right[i]+1,this_genes$Left[(i+1)]-1, this_genes$Left[(i+1)]-this_genes$Right[i]-1),index=k)
      }
      k <- k+1
    }
    # Insert a final "intergenic" region if the final coordinate of the chromosome does not contain a gene
    if(!this_genes$Right[nrow(this_genes)]==chrlength){
      new_this_genes <- insertRow(new_this_genes,newrow=c("IG","",1,this_genes$Right[nrow(this_genes)]+1, chrlength, chrlength-this_genes$Right[nrow(this_genes)]),index=k+1)
    }
    rownames(new_this_genes) <- 1:nrow(new_this_genes)
  }
  new_this_genes[,3:6] <- sapply(new_this_genes[,3:6],function(x) (as.numeric(x)))
  return(new_this_genes)
}
