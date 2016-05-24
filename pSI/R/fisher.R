#'@title 
#'Fisher's Exact Test (Single Cell Type & pSI Threshold)

fisher<-function(pSI, candidate.genes, thresh, cell, total){
  #creating index of genes and column binding it onto the existing pSI data frame
  pSI<-cbind(c(1:length(pSI[,1])),pSI) 
  
  #compensating for cell/sample column number because of the column bind just performed
  cell <- cell+1
  
  #this identifies, for a given cell type (column) from the pSI table, the genes which have pSI values less than a certain threshold
  hits <- pSI[which(pSI[,cell]<thresh),1]
  
  #as long as there is one hit, calculate number of overlapping genes below 
  if(!is.na(hits[1])){ # 7_26_13 
    #hits is the list of row number for an gene with pSI<thresh. It runs through all these hits...
    for (i in 1:length(hits)){   
      k<-hits[i]
      #and asks if they are listed as a hit on your gene list, and makes a list of the overlap called outy.
      if (candidate.genes[candidate.genes[,1]==k, 2]  >0){
        if (exists("outy")){
          outy<-rbind(outy, candidate.genes[  candidate.genes[,1]==k, ]  )
        } else {                 
          outy     <- candidate.genes[  candidate.genes[,1]==k, ]
        }
      }
    }
  }    
  
  if(exists("outy")){
    if(length(outy)==2){
      a<-1
    }else{
      #number overlapping is equal to this list
      a<-length(outy[,1]) 
    }                     
    
  } else {
    a<-0
  } 
  
  #number of human candidate genes not overlapping
  b<-sum(candidate.genes[,2]>0) - a    
  #number of cell type spec genes not overlapping
  c<-length(hits) -a  
  #number of other genes.
  d<-total -a -b -c
  
  #format expected by fisher.test function.
  counts<-matrix(c(a,c,b,d),nrow=2)
  #this is the p-value, which is all you get otherwise. 
  fisher.test(counts, alternative="greater")[1]
  
}