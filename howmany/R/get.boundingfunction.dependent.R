"get.boundingfunction.dependent" <-
function(X,Y,alpha,test,alternative,at=(1:1000)/1000,n.permutation=round(20/alpha))
  {
    
    m <- ncol(X)
    n <- nrow(X) 
     
    cat(" ... ", n.permutation," permutations\n ")
    Testall <- matrix(rep(0,m*n.permutation),nrow=n.permutation)
    Quantile <- matrix(rep(0,m*n.permutation),nrow=n.permutation)

    for (perm in 1:n.permutation)
      { 
        if(round(perm/50)==perm/50) cat( perm, "  ")
        Ynew <- Y[sample(1:n,n)]
        for (p in 1:m){
          Testall[perm,p] <- test(X[Ynew==0,p],X[Ynew==1,p],alternative=alternative)$p.value
        }
      }
    cat("  \n")
    Quantile <- Testall
    #for (p in 1:m) Quantile[,p] <- sample(Quantile[,p],n.permutation)
    for (perm in 1:n.permutation) Quantile[perm,] <- sort(Quantile[perm,])
    for (p in 1:m) Quantile[,p] <- sort(Quantile[,p])
    
    quanbeta <- 1
    too.high <- 1

  
    while(too.high)
      {
        quanbeta <- quanbeta+1
        
        boundingjump <- rep(0,m)
        boundingjump <- Quantile[quanbeta,]
        
        count <- 0
        for (perm in 1:n.permutation)
          {
            count <- count  +  min(1,   sum(Testall[perm,]<boundingjump)  )
          }
        
        if(count>alpha*n.permutation)
          {
            too.high <- 0
            quanbeta <- quanbeta-1
            boundingjump <- Quantile[quanbeta,]
          }
      }

    boundingfunction <- numeric(length(at))
    for (i in 1:length(at))
      {
        boundingfunction[i] <- sum(boundingjump<=at[i])
      }
       
    return(boundingfunction)
        
  }

