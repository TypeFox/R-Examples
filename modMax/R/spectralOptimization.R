#execute general spectral optimization algorithm with optional refinement
spectralOptimization <- function(adjacency, numRandom=0,initial=c("general","own"),
                                 refine=FALSE){
  
  initial <- match.arg(initial)
  #refine <- match.arg(refine)
  
  if(initial=="own"){
    C_initial <- adjacency[,length(adjacency[1,])]
    adjacency <- adjacency[,-length(adjacency[1,])]
  }
  else{
    C_initial <- seq(0,0,length.out=length(adjacency[,1]))
  }
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  res <- callSpectralOptimization(network,initialC=C_initial,refine=refine)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callSpectralOptimization(network,initialC=seq(0,0,length.out=vcount(network)),
                                            refine=refine)
      rndResult[vcount(rndNetwork)+1] <- round(rndResult[vcount(rndNetwork)+1],2)
      randomResults <- rbind(randomResults,rndResult)
    }
    result <- generateOutput(res,randomResults,random=TRUE)
  }
  else{
    result <- generateOutput(res)
  }
  
  return(result)
}

#execute multi way division approach according to Wang et al.
multiWay <- function(adjacency, numRandom=0, maxComm=length(adjacency[1,])){
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  res <- callMultiWay(network,maxComm)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callMultiWay(rndNetwork,maxComm)
      rndResult[vcount(rndNetwork)+1] <- round(rndResult[vcount(rndNetwork)+1],2)
      randomResults <- rbind(randomResults,rndResult[vcount(rndNetwork)+1])
    }
    result <- generateOutput(res,randomResults,random=TRUE)
  }
  else{
    result <- generateOutput(res)
  }
  
  return(result)
}

#execute spectral1
spectral1 <- function(adjacency, numRandom=0, maxComm=(length(adjacency[1,])-1)){
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  res <- callSpectral1(network,maxComm)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callSpectral1(rndNetwork,maxComm)
      rndResult[vcount(rndNetwork)+1] <- round(rndResult[vcount(rndNetwork)+1],2)
      randomResults <- rbind(randomResults,rndResult[vcount(rndNetwork)+1])
    }
    result <- generateOutput(res,randomResults,random=TRUE)
  }
  else{
    result <- generateOutput(res)
  }
  
  return(result)
}

#execute spectral-2 algorithm
spectral2 <- function(adjacency, numRandom=0, maxComm=(length(adjacency[1,])-1)){
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  res <- callSpectral2(network,maxComm)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callSpectral2(rndNetwork,maxComm)
      randomResults <- rbind(randomResults,rndResult)
    }
    result <- generateOutput(res,randomResults,random=TRUE)
  }
  else{
    result <- generateOutput(res)
  }
  
  return(result)
}

####functions for execution of the different spectral optimization algorithms
#general spectral optimization with optional refinement
callSpectralOptimization <- function(network, initialC=seq(0,0,length.out=vcount(network)),
                                 refine){
  
  if(initialC[1]==0){
    C <- initialWholeGraph(network)
    Q <- 0
  }
  else{
    C <- initialC
    Q <- calculateQ(network,C)
  }
 
  modMatrix <- calculateModMatrix(network)
  
  changed <- 1
  while(changed==1){
    changed <- 0
    C_new <- C
    for(c in 1:max(C)){
      #print(c)
      
      sizeComm <- length(which(C==c))
      indices <- which(C==c)
      B_c <- matrix(,ncol=sizeComm,nrow=sizeComm)
      for(i in 1:sizeComm){
       
        currI <- indices[i]
        for(j in i:sizeComm){
          currJ <- indices[j]
          B_c[i,j] <- modMatrix[currI,currJ]
          if(C[currI]==C[currJ]){
            vertices <- which(C==C[currI])
            substract <- sum(modMatrix[currI,vertices])
            B_c[i,j] <- B_c[i,j]-substract
          }
          B_c[j,i] <- B_c[i,j]
        }  
      }

      eigenvectors <- eigen(B_c)
      value <- max(eigenvectors$values)
      index <- which(eigenvectors$values==value)
      vector <- eigenvectors$vectors[,index]
      
      if(value>0&!is.complex(value)&(length(which(vector>0))>0)&(length(which(vector<0))>0)){
        s <- seq(0,0,length.out=sizeComm)
        for(i in 1:length(s)){
          if(vector[i]<0){
            s[i] <- -1
          }
          else{
            s[i] <- 1
          }
        }
        
        deltaQ <- s%*%B_c%*%s/(2*sum(get.adjacency(network,attr="weight")))

        if(deltaQ > 0){
          C_tmp <- seq(1,1,length.out=length(s))
          for(i in 1:length(s)){
            if(s[i]==(-1)){
              C_tmp[i]<-2
            }
          }
          changed <- 1
          
          if(refine==TRUE){
            subgraph <- induced.subgraph(network,indices)
            Q_tmp <- calculateQ(subgraph,C_tmp)
            result <- kernighanLinRefinement(subgraph,C_tmp,Q_tmp)
            C_tmp <- result[1:vcount(subgraph)]
            if(length(which(C_tmp==1))==0|length(which(C_tmp==2))==0){
              changed <- 0
            }
          }
          currMax <- max(C_new)
          for(i in 1:length(C_tmp)){
            currI <- indices[i]
            if(C_tmp[i]==2){
              C_new[currI]=currMax+1
            }
          }
          C_new <- updateC(C_new)
        }
      }
    }
    C <- C_new
  }
  
  Q <- calculateQ(network,C)
  result <- c(C,Q)
  
  return(result)
}

#multiway division approach according to Wang et al.
callMultiWay <- function(network,maxComm){
  
  C_best <- initialWholeGraph(network)
  Q_best <- 0
  
  modMatrix <- calculateModMatrix(network)
  eigens <- eigen(modMatrix)
  eigenvalues <- eigens$values
  eigenvalues <- sort(eigenvalues, decreasing=TRUE)
  indices <- NULL
  for(e in eigenvalues){
    indices <- c(indices,which(eigens$values==e))
  }
  eigenvectors <- eigens$vectors[,indices]
  
  posEigen <- length(which(eigenvalues>0))
  if(posEigen<(maxComm-1)){
    maxComm <- posEigen+1
  }
  
  eigenvalues <- eigenvalues[1:maxComm]
  eigenvectors <- eigenvectors[,1:maxComm]
  
  for(h in 2:maxComm){
    #percentage <- (h-1)/maxComm*100
    #print(paste(percentage,"%",sep=""))
    
    U_h <- eigenvectors[,1:(h-1)]
    val_h <- eigenvalues[1:(h-1)]

    #compute V
    V <- NULL
    if(h>2){
      for(i in 1:length(U_h[,1])){
        row <- NULL
        for(k in 1:length(U_h[1,])){
          currEl <- sqrt(val_h[k])*U_h[i,k]
          row <- cbind(row,currEl)
        }
        V <- rbind(V,row)
      }
    }
    else{
      V <- sapply(U_h, function(x) x=sqrt(val_h)*x)
    }
    
    S <- NULL
    for(i in 1:vcount(network)){
      s_i <- seq(1,1,length.out=(h-1))
      if(h>2){
        for(k in 1:(h-1)){
          if(V[i,k]<0){
            s_i[k] <- -1
          }
        }
      }
      else{
        if(V[i]<0){
          s_i <- -1
        }
      }
      
      S <- rbind(S,s_i)
    }

    groups <- permutations(n=2,r=h-1,v=c(-1,1),repeats.allowed=TRUE)
    
    Ygroups <- seq(0,0,length.out=vcount(network))
    numGroups <- seq(0,0,length.out=length(groups[,1]))
    for(i in 1:length(groups[,1])){
      indices <- NULL
      if(h>2){ 
        for(j in 1:length(S[,1])){
          if(isTRUE(all.equal(S[j,],groups[i,]))){
            indices <- c(indices,j)
          }
        }
      }
      else{
        indices <- which(S==groups[i])
      }

      Ygroups[indices] <- i
      numGroups[i] <- length(indices)
    }

    n <- length(numGroups)
    
    part <- n-h+1
    hGroups <- sort(numGroups)[part:n]
    maxGroups <- NULL
    for(x in hGroups){
      maxG <- which(numGroups==x)
      maxGroups <- c(maxGroups,maxG)
    }
    maxGroups <- unique(maxGroups)
    if(length(maxGroups)>h){
      diff <- length(maxGroups)-h
      maxGroups <- maxGroups[-diff]
    }
        
    remainingVertices <- seq(1,vcount(network),1)
    C_h <- seq(0,0,length.out=vcount(network))
    Y <- NULL

    for(j in maxGroups){
      Y_j <- 0
      for(i in 1:vcount(network)){
        if(Ygroups[i]==j){
          if(h>2){
            Y_j <- Y_j+V[i,]
          } else{
            Y_j <- Y_j+V[i]
          }
          remainingVertices <- remainingVertices[remainingVertices !=i]
          C_h[i]<-j
        }  
      }
      Y <- cbind(Y,Y_j)
    }
    
    for(m in remainingVertices){
      Yy_max <- -Inf
      j_max <- -1
      for(j in 1:h){
        if(h>2){
          tmp <- Y[,j]%*%V[m,]
        }
        else{
          tmp <- Y[1,j]*V[m]
        }
        if(tmp>Yy_max){
          Yy_max <- tmp
          j_max <- j
        }     
      }
      C_h[m] <- maxGroups[j_max]
    }

    C_h <- updateC(C_h)
    Q_h <- calculateQ(network,C_h)
    
    if(Q_h>Q_best){
      Q_best <- Q_h
      C_best <- C_h
    }
  }
  
  result <- kernighanLinRefinement(network,C_best,Q_best)
  return(result)
}

#spectral-1 algorithm
callSpectral1 <- function(network,maxComm){
  
  U_K <- computeUk(network,maxComm)
  maxComm <- length(U_K[1,])
  
  C <- initialWholeGraph(network)
  Q <- 0
  for(k in 2:maxComm){
    u_k <- U_K[,1:k]
    #print(dim(u_k))
    for(j in 1: length(u_k[,1])){
      squareVec <- sapply(u_k[j,], function(x) x=x^2)
      normvec <- sqrt(sum(squareVec))
      u_k[j,] <- sapply(u_k[j,], function(x) x=x/normvec)
    }
    #print(k)
    C_k <- kmeans(u_k,k)$cluster
    Q_k <- calculateQ(network,C_k)
    
    if(Q_k>Q){
      Q <- Q_k
      C <- C_k
    }
  }
  
  result <- c(C,Q)
  return(result)
}

#spectral-2 algorithm
callSpectral2 <- function(network,maxComm){
  
  U_K <- computeUk(network,maxComm)
  maxComm <- length(U_K[1,])
  C <- initialWholeGraph(network)
  Q <- 0
  split <- 1
  exists <- 0
  k <- 2
  while(k<=maxComm & split==1){
    C_new <- C
    Q_new <- Q
    split <- 0
    for(i in 1:max(C)){
      #print(C)
      vertices <- which(C==i)
      
      if(exists==0){
        u_k <- U_K[,1:k]
        for(j in 1: length(u_k[,1])){
          squareVec <- sapply(u_k[j,], function(x) x=x^2)
          normvec <- sqrt(sum(squareVec))
          u_k[j,] <- sapply(u_k[j,], function(x) x=x/normvec)
        }
        exists <- 1
      }
      if(length(vertices)>1){
        #print(vertices)
        u_kc <- u_k[vertices,]
        #print(u_kc)
        if(length(vertices)>2){
          tmpComm <- kmeans(u_kc,2)$cluster
        }
        else{
          tmpComm <- c(1,2)
        }

        C_tmp <- C_new
        for(j in 1:length(tmpComm)){
          vertex <- vertices[j]
          if(tmpComm[j]==2){
            C_tmp[vertex]<-max(C_new)+1
          }
        }
        Q_tmp <- calculateQ(network,C_tmp)
        if(Q_tmp>Q_new){
          Q_new <- Q_tmp
          C_new <- C_tmp
          k <- k+1
          split <- 1
          exists <- 0
        }
      }
    }
  
    C <- C_new
    Q <- Q_new
  }
  
  result <- c(C,Q)
  return(result)
}

computeUk <- function(network,maxComm){
  adjacency <- get.adjacency(network,attr="weight")
  D <- matrix(0,nrow=vcount(network),ncol=vcount(network))
  for(i in 1:vcount(network)){
    D[i,i] <- sum(adjacency[i,])
  }
  M <- solve(D)%*%adjacency
  eigens <- eigen(M)
  eigenvalues <- eigens$values
  eigenvalues <- sort(eigenvalues, decreasing=TRUE)
  indices <- NULL
  for(e in eigenvalues){
    indices <- cbind(indices,which(eigens$values==e))
  }
  eigenvectors <- eigens$vectors[,indices]
  
  eigenvalues <- eigenvalues[1:maxComm]
  eigenvalues <- eigenvalues[Re(eigenvalues)>0]
  dimension <- length(eigenvalues)
  eigenvectors <- eigenvectors[,1:dimension]
  
  return(eigenvectors)
}