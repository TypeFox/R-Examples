
populate <- function(input,dim,simslist=FALSE,samples=NULL){
  
  if(!simslist){
    
    charinds <- sub(".*\\[(.*)\\].*", "\\1", names(input), perl=TRUE) 
  
    fill <- array(NA,dim=dim)
    
    for (i in 1:length(input)){
      
      ind <- lapply(strsplit(charinds[i], ','), as.integer)[[1]]
      fill[matrix(ind,1)] <- input[i]
    
    }
  } else {
    
    charinds <- sub(".*\\[(.*)\\].*", "\\1", colnames(input), perl=TRUE) 
    
    fill <- array(NA,dim=c(samples,dim))
    
    for (i in 1:length(charinds)){
      
      #ind <- lapply(strsplit(charinds[i], ','), as.integer)[[1]]
      
      eval(parse(text=paste('fill[','1:',samples,',',charinds[i],']','<- input[,i]',sep="")))
      
    }
  }

  return(fill)
  
}



