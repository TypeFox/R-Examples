patch_age <- function(span)
{
  s1 <- list()
  for(i in 1:length(span)){
    
    s1[[i]] <- span[[i]][,8]
  }
  
  s2 <- unique(unlist(s1))
  max1 <- length(s2)
     
  mat1 <- as.data.frame(matrix(nrow=max1,ncol=length(span)))
  rownames(mat1) <- s2
  
  for(y in 1:ncol(mat1)){
    
    s0 <- as.vector(s1[[y]])
    l_s0 <- length(s0)
    
    if(l_s0!=0){
    for (v in 1:l_s0){
      val1 <- s0[v]
      if (val1%in% rownames(mat1)) mat1[rownames(mat1)==val1,y] <- 1 
    }}
  } 
  
  mat1[is.na(mat1)] <- 0 
  
  for(b in 1:nrow(mat1)){
    
    l2 <- as.numeric(mat1[b,])
    #Cumulative sum
    l3 <- cumsum(l2)
    l3 <- l3*l2
    
    mat1[b,] <- l3
  }
  
  for(m in 1:nrow(mat1)){
    
    l1 <- as.numeric(mat1[m,])
    
    for(n in 1:length(l1)){
      
      if(l1[n]!=0)break
      
      else l1[n]<-NA
    } 
    
    mat1[m,] <- l1	
    
  }
  
  span1 <- list()
  
  for(d in 1:length(span)){
    
    r1 <- span[[d]]
    
    col1 <- mat1[,d]
    
    col2 <- col1[!is.na(col1)]
    
    col3 <- col2[col2 != 0 ] 
    
    r2 <- cbind(r1,col3)
    
    colnames(r2)[9] <- "age"
    
    span1[[d]] <- r2
    
  }  
  
  return(span1)
}