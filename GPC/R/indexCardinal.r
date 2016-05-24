indexCardinal <- function(d,p,m,index){
  if (d == 1){
    index <- matrix(seq(0,p),nrow = 1,ncol = p+1)      
  } else {
    if (nargs() == 2){
      m = getM(d,p)
      index = matrix(0,nrow=d,ncol=m)
    }
    if (m <= d+1){
      for (i in 1:d){index[i,i+1] <- 1.0}
    } else {
      m2vec <- vector(length=p+1)
      for (ii in 0:p){
        m2vec[ii+1] <- getM(d,ii)} # getM(ii,d)
      cardp <- sum(1*(m>m2vec))
      cardp -> pm
      mnn <- (m - m2vec[pm])
      for (nn in d:2){
        m1vec <- m2vec                                
        for (ii in 0:(pm+1)){
          m2vec[ii+1] <- getM(nn-1,ii)} # getM(ii,nn-1)  		
        pnn <- sum(1*(mnn > m2vec))
        index[d-nn+1,m] <- (pm - pnn) 
        mnn <- 1*(mnn - m2vec[pnn])
        pm <- pnn
      }
      index[d,m]=cardp-sum(index[,m])
      index = indexCardinal(d,p,m-1,index)	
    }    
  }
  return(index)	
}