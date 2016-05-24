#execute genetic algorithm
geneticAlgorithm <- function(adjacency, numRandom=0,initial=c("general","cluster","own"),p,g,mutRat=0.5,
                             crossOver=0.2,beta=0.1,alpha=0.4,n_l=4,local=FALSE){
  
  initial <- match.arg(initial)
  
  if(initial=="own"){
    C_initial <- adjacency[,length(adjacency[1,])]
    adjacency <- adjacency[,-length(adjacency[1,])]
  }
  else{
    C_initial <- seq(0,0,length.out=length(adjacency[,1]))
  }
  
  network <- graph.adjacency(adjacency,mode="undirected",weighted=TRUE)
  
  res <- callGeneticAlgorithm(network,initialC=C_initial,initial=initial,
                                 p=p,g=g,mutRat=mutRat,crossOver=crossOver,
                                 beta=beta,alpha=alpha,n_l=n_l,local=local)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callGeneticAlgorithm(rndNetwork,
                                        initialC=seq(0,0,length.out=vcount(rndNetwork)),initial=initial,
                                        p=p,g=g,mutRat=mutRat,crossOver=crossOver,
                                        beta=beta,alpha=alpha,n_l=n_l,local=local)
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

####functions for execution of genetic algorithm
callGeneticAlgorithm <- function(network,initialC=seq(0,0,length.out=vcount(network)),
                                                      initial,p,g,mutRat,crossOver,beta,alpha,n_l,local){
  
  P <- NULL
  for(i in 1:p){
    if(initial=="general"){
      C <- initialSingletons(vcount(network))
    }
    else if(initial=="cluster"){
      C <- initialSingletons(vcount(network))
      iter <- round(alpha*vcount(network),0)
      for(init in 1:iter){
        vertex <- sample(1:vcount(network),1)
        neighs <- neighborhood(network,1,vertex)
        neighbors <- unlist(neighs)
        for(v in neighbors){
          C[v]=C[vertex]
        }
      }
      C <- updateC(C)
    }
    else if(initialC[1]!=0){
      C <- initialC
    }

    P <- rbind(P,C)
  }
  for(i in 1:g){
    #percentage <- (i-1)/g*100
    #print(paste(percentage,"%",sep=""))
    
    #apply fitness function to the chromsomes
    fitness <- calculateFitness(network,P)
    
    fitness <- fitness[order(fitness[,1],decreasing=TRUE),]
    pop <- fitness[1:p,2]
    newP <- P[pop,]
    fract <- round(beta*p,0)
    save <- fitness[1:fract,2]
    saveP <- P[save,]
    
    #apply cross over operation
    indices <- seq(1,length(newP[,1]),2)
    for(k in indices){
      if(k<length(newP[,1])){
        src <- sample(c(k,k+1),1)
        if(src == k){
          dest <- k+1
        } else{
          dest <- k
        }
        crossOverNum <- round(crossOver*vcount(network),0)
        for(j in 1:crossOverNum){
          vertex <- sample(1:vcount(network),1)
          cluster <- newP[src,vertex]
          vertices <- which(newP[src,]==cluster)
          for(v in vertices){
            newP[dest,v]=cluster
          }
        }
      }
    }
    
    #apply point mutation
    mutNum <- round(mutRat*vcount(network),0)
    for(h in 1:mutNum){
      chromo <- sample(1:length(newP[,1]),1)
      vertices <- sample(1:vcount(network),2,replace=FALSE)
      v1 <- vertices[1]
      v2 <- vertices[2]
      newP[chromo,v2] <- newP[chromo,v1]
    }
    for(c in 1:length(newP[,1])){
      newP[c,] <- updateC(newP[c,])
    }
    
    #apply local search operator
    if(local==T){
      newP <- localSearchOperator(network,P,n_l)
    }
    P <- rbind(newP,saveP)
  }
  fitness <- calculateFitness(network,P)
  index <- which.max(fitness[,1])
  C_best <- P[index,]
  Q_best <- max(fitness[,1])
  result <- c(C_best,Q_best)
  
  return(result)
}

#calculate the fitness values for a population
calculateFitness <- function(network,P){
  fitness <- NULL
  for(k in 1:length(P[,1])){
    Q <- calculateQ(network,P[k,])
    currQ <- cbind(Q,k)
    fitness <- rbind(fitness,currQ)
  }
  return(fitness)
}

#apply local search operator
localSearchOperator <- function(network,P,n_l){
  
  A <- get.adjacency(network,attr="weight")
  newP <- NULL
  for(k in 1:length(P[,1])){
    chromo <- P[k,]
    newChromo <- chromo
    bestChromo <- chromo
    fitness <- calculateQ(network,chromo)
    bestFitness <- fitness
    for(i in 1:n_l){
      row <- sample(1:length(A[,1]),1)
      for(q in 1:length(A[row,])){
        if(A[row,q]!=0){
          newChromo[q]<-chromo[row]
        }
      }
      newFitness <- calculateQ(network,newChromo)
      if(newFitness>bestFitness){
        bestChromo <- newChromo
        bestFitness <- newFitness
      }
    }
    if(bestFitness>fitness){
      bestChromo <- updateC(bestChromo)
      newP <- rbind(newP,bestChromo)
    } else{
      newP <- rbind(newP,chromo)
    }
  }
  
  return(newP)
}