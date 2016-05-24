#simulated annealing algorithm according to Medus/Massen&Doye
simulatedAnnealing <- function(adjacency, numRandom=0, 
                               initial=c("general","random","greedy","own"),beta=length(adjacency[1,])/2,alpha=1.005,fixed){
  
  initial <- match.arg(initial)
  
  if(initial=="own"){
    C_initial <- adjacency[,length(adjacency[1,])]
    adjacency <- adjacency[,-length(adjacency[1,])]
  }
  else{
    C_initial <- seq(0,0,length.out=length(adjacency[,1]))
  }
  
  network <- graph.adjacency(adjacency,mode="undirected",weighted=TRUE)  
  res <- callSimulatedAnnealing(network,initialC=C_initial,initial=initial,beta=beta,alpha=alpha,fixed=fixed)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callSimulatedAnnealing(rndNetwork,
                                          initialC=seq(0,0,length.out=vcount(rndNetwork)),
                                          initial=initial,beta=beta,alpha=alpha,fixed=fixed)
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

#simulated annealing algorithm according to Guimera & Amaral
saIndividualCollectiveMoves <- function(adjacency,numRandom=0,initial=c("general","own"),beta=length(adjacency[1,])/2,alpha=1.005,fixed=25,numIter=1.0){
  
  initial <- match.arg(initial)
  
  if(initial=="own"){
    C_initial <- adjacency[,length(adjacency[1,])]
    adjacency <- adjacency[,-length(adjacency[1,])]
  }
  else{
    C_initial <- seq(0,0,length.out=length(adjacency[,1]))
  }
  
  network <- graph.adjacency(adjacency,mode="undirected",weighted=TRUE)  
  res <- callSaIndividualCollectiveMoves(network,initialC=C_initial,beta=beta,alpha=alpha,fixed=fixed,numIter=numIter)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callSaIndividualCollectiveMoves(rndNetwork,
                                                   initialC=seq(0,0,length.out=vcount(rndNetwork)),
                                                   beta=beta,alpha=alpha,fixed=fixed,numIter=numIter)
      result[vcount(network)+1] <- round(result[vcount(network)+1],2)
      randomResults <- rbind(randomResults,rndResult[vcount(rndNetwork)+1])
    }
    result <- generateOutput(res,randomResults,random=TRUE)
  }
  else{
    result <- generateOutput(res)
  }
  
  return(result)
}

####functions called when executing of simulated annealing algorithms
callSimulatedAnnealing <- function(network, 
                                   initialC=seq(0,0,length.out=vcount(network)),
                                    initial=c("general","random","greedy","own"),
                                    beta, alpha, fixed){
  initial <- match.arg(initial)
  if(initial=="general"){
    C <- initialSingletons(vcount(network))
  }
  else if(initial=="random"){
    C <- initialRandom(network)
  }
  else if(initial=="greedy"){
    C <- initialGreedy(network)
  }
  else if(initialC[1]!=0){
    C <- initialC
  }
  
  Q <- calculateQ(network,C)
  C_best <- C
  Q_best <- Q
  
  steps <- 0
  while(steps<fixed){
    #percentage <- steps/fixed*100
    #print(paste(percentage,"%",sep=""))
    
    tmpResult <- individualNodeMovement(network,C,Q,beta)
    Q_new <- tmpResult[vcount(network)+1]
    C_new <- tmpResult[1:vcount(network)]
    accepted <- tmpResult[vcount(network)+2]

    if(Q_new > Q_best){
      Q_best <- Q_new
      C_best <- C_new
    }
    C <- C_new
    Q <- Q_new
    if(accepted==0){
      steps <- steps+1
    } else{
      steps <- 0
    }
    beta <- beta*alpha
  }
  
  result <- c(C_best, Q_best)
  return(result)
}

#simulated annealing algorithm according to Guimera & Amaral
callSaIndividualCollectiveMoves <- function(network,initialC=seq(0,0,length.out=vcount(network)),
                                            beta=vcount(network)/2,alpha=1.005,fixed,numIter){
  
  if(initialC[1]==0){
    C <- initialGreedy(network)
  }
  else{
    C <- initialC
  }

  Q <- calculateQ(network,C)
  beta_initial <- beta
  
  C_best <- C
  Q_best <- Q
  
  n <- vcount(network) 
  steps <- 0

  while(steps < fixed){
    #percentage <- steps/fixed*100
    #(paste(percentage,"%",sep=""))
    #print(C_best)
    
    for(a in 1:numIter*n^2){
      result <- individualNodeMovement(network, C, Q, beta)
      Q_new <- result[vcount(network)+1]
      C_new <- result[1:vcount(network)]
      accepted <- result[vcount(network)+2]
      if(Q_new > Q_best){
        Q_best <- Q_new
        C_best <- C_new
      }
      Q <- Q_new
      C <- C_new
      if(accepted==0){
        steps <- steps+1
      } else{
        steps <- 0
      }
    }
    for(b in 1:numIter*n){
      if(max(C)>1){
        split <- sample(c(0,1),1)
      }
      else{
        split <- 1
      }
      
      if(split==0){
        result <- mergeCommunities(network,C,Q,beta)
      }
      else{
        result <- splitCommunities(network,C,Q,beta_initial,alpha,beta)
      }
      Q_new <- result[vcount(network)+1]
      C_new <- result[1:vcount(network)]
      accepted <- result[vcount(network)+2]
      if(Q_new > Q_best){
        Q_best <- Q_new
        C_best <- C_new
      }
      Q <- Q_new
      C <- C_new
      if(accepted==0){
        steps <- steps+1
      } else{
        steps <- 0
      }
    }
    beta = beta*alpha
  }

  result <- c(C_best, Q_best)
  return(result) 
}