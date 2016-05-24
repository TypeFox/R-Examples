clusterscore <-
  function(Z)
  {
    #Z is a centered and scaled numerical matrix
    Z<-as.matrix(Z)
    n<-nrow(Z)
    Ztilde <-Z/sqrt(n)
    e <- svd(Ztilde)
    f <- e$u[,1]*e$d[1]*sqrt(n) #first principal component
    sv<-e$d[1]	#standard deviation of f i.e. first singular value
    v <- e$v[,1]
    return(list(f=f,sv=sv,v=v))	
  }

