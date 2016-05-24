#execute the greedy algorithm
greedy <- function(adjacency,numRandom=0, 
                   q=c("general","danon","wakita1","wakita2","wakita3"),
                   initial=c("general","prior","walkers","subgraph","adclust","own"),
                   randomized=0,
                   refine=c("none","complete","fast","kernighan"),
                   coarse=0){
  
  q <- match.arg(q)
  initial <- match.arg(initial)
  refine <- match.arg(refine)
  
  if(initial=="own"){
    C_initial <- adjacency[,length(adjacency[1,])]
    adjacency <- adjacency[,-length(adjacency[1,])]
  }
  else{
    C_initial <- seq(0,0,length.out=length(adjacency[,1]))
  }
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)  
  res <- callGreedy(network,initialC=C_initial,q=q,initial=initial,
                    randomized=randomized, refine=refine,coarse=coarse)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callGreedy(rndNetwork,initialC=seq(0,0,length.out=vcount(rndNetwork)),q=q,initial=initial,
                              randomized=randomized, refine=refine,coarse=coarse)
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

#execute the RG+ approach
rgplus <- function(adjacency,numRandom=0,z,randomized){
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  res <- callRgplus(network,z,randomized)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom !=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callRgplus(rndNetwork,z=z,
                              randomized=randomized)
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

#execute the msg-vm approach
msgvm <- function(adjacency,numRandom=0,initial=c("general","own"), parL){
  
  initial <- match.arg(initial)
  
  if(initial=="own"){
    C_initial <- adjacency[,length(adjacency[1,])]
    adjacency <- adjacency[,-length(adjacency[1,])]
  }
  else{
    C_initial <- seq(0,0,length.out=length(adjacency[,1]))
  }
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  res <- callMsgvm(network,initialC=C_initial,parL=parL)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom !=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callMsgvm(rndNetwork,initialC=seq(0,0,length.out=vcount(network)),parL=parL)
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

#execute the cd algorithm
cd <- function(adjacency, numRandom=0,initial=c("general","own"),maxC=length(adjacency[,1]),iter,p){
  
  initial <- match.arg(initial)
  
  if(initial=="own"){
    C_initial <- adjacency[,length(adjacency[1,])]
    adjacency <- adjacency[,-length(adjacency[1,])]
  }
  else{
    C_initial <- seq(0,0,length.out=length(adjacency[,1]))
  }
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  res <- callCD(network,initialC=C_initial,maxC=maxC,iter=iter,p=p)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callCD(rndNetwork,initialC=seq(0,0,length.out=vcount(network)),maxC=maxC,iter=iter,p=p)
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

#execute the louvain algorithm
louvain <- function(adjacency, numRandom=0, initial=c("general","own")){
  
  initial <- match.arg(initial)
  
  if(initial=="own"){
    C_initial <- adjacency[,length(adjacency[1,])]
    adjacency <- adjacency[,-length(adjacency[1,])]
  }
  else{
    C_initial <- seq(0,0,length.out=length(adjacency[,1]))
  }
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  res <- callLouvain(network,initialC=C_initial)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callLouvain(rndNetwork,initialC=seq(0,0,length.out=vcount(rndNetwork)))
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

#execute greedy approach based on vertex similarity and using hybrid merge
vertexSim <- function(adjacency, numRandom=0, frac=0.5){
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  res <- callVertexSim(network,frac=frac)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callVertexSim(rndNetwork,frac=frac)
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

#execute mome algorithm
mome <- function(adjacency, numRandom=0){
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  
  res <- callMome(network)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callMome(rndNetwork)
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

#####Greedy functions#####
#Greedy algorithm and its modifications
callGreedy <- function(network,initialC=seq(0,0,length.out=vcount(network)),
                       q=c("general","danon","wakita1","wakita2","wakita3"),
                       initial=c("general","prior","walkers","subgraph","adclust"),
                       randomized=0,
                       refine=c("none","complete","fast","kernighan","adclust"),
                       coarse=0){
  
  #initialization
  q <- match.arg(q)
  initial <- match.arg(initial)
  refine <- match.arg(refine)
  
  if(initialC[1]!=0){
    C <- initialC
    a <- calculateA(network,C)
    DeltaQ <- calculateDeltaQ(network,C)
  }
  else if(initial=="general"){
    C <- initialSingletons(vcount(network))
    DeltaQ <- matrix(,nrow=max(C),ncol=max(C))
    a <- seq(0,0,length.out=max(C))
    for(i in 1:length(a)){
      a[i] <- getDegree(network,i)/(sum(get.adjacency(network,attr="weight")))  
    }
    for(i in 1:(max(C)-1)){
      for(j in (i+1):max(C)){
        if(are.connected(network,i,j)){
          DeltaQ[i,j] <- 2/(sum(get.adjacency(network,attr="weight")))-
            2*getDegree(network,i)*getDegree(network,j)/(sum(get.adjacency(network,attr="weight"))*sum(get.adjacency(network,attr="weight")))
          DeltaQ[j,i] <- DeltaQ[i,j]
        }
      }  
    }
  }
  else if(initial=="prior"){
    C <- initialPriorKnowledge(network)
    a <- calculateA(network,C)
    DeltaQ <- calculateDeltaQ(network,C)
  }
  else if(initial=="walkers"){
    C <- initialRandomWalkers(network)
    a <- calculateA(network,C)
    DeltaQ <- calculateDeltaQ(network,C)
  }
  else if(initial=="subgraph"){
    C <- initialSubgraphSim(network)
    a <- calculateA(network,C)
    DeltaQ <- calculateDeltaQ(network,C)
  }
  else if(initial=="adclust"){
    
    refine <- "adclust"
    coarse <- 0
    
    C <- initialSingletons(vcount(network))
    Q <- calculateQ(network,C)
    tmp <- fastGreedyRefinement(network,C,Q)
    C <- tmp[1:vcount(network)]
    
    a <- calculateA(network,C)
    DeltaQ <- calculateDeltaQ(network,C)
  }
  
  if(coarse !=0){
    numComm <- max(C)
  }
  
  Q <- calculateQ(network,C)
  Q_best <- Q
  C_best <- C
  
  while(max(C)>1&length(DeltaQ)>1&length(which(!is.na(DeltaQ)))>0){
    #percentage <- (vcount(network)-max(C))/vcount(network)*100
    #print(paste(percentage,"%",sep=""))
    
    if(q=="general"){
      maxEl <- which(DeltaQ==max(DeltaQ,na.rm=TRUE),arr.ind=TRUE)[1,]
    }
    else if(q=="danon"){
      DeltaQDanon <- DeltaQ
      for(i in 1:length(DeltaQ[,1])){
        DeltaQDanon[i,] <- sapply(DeltaQDanon[i,], function(x) x=x/a[i])
      }
      maxEl <- which(DeltaQDanon==max(DeltaQDanon,na.rm=TRUE),arr.ind=TRUE)[1,]
    }
    else if(q=="wakita1"){
      DeltaQWakita <- DeltaQ
      for(i in 1:(length(DeltaQ[,1])-1)){
        sizeCommI <- length(na.omit(DeltaQ[i,]))
        for(j in (i+1):length(DeltaQ[,1])){
          if(!is.na(DeltaQ[i,j])){
            sizeCommJ <- length(na.omit(DeltaQ[j,]))
            consol <- min(sizeCommI/sizeCommJ,sizeCommJ/sizeCommI)
            DeltaQWakita[i,j] <- DeltaQ[i,j]*consol
          }
        }  
      }
      maxEl <- which(DeltaQWakita==max(DeltaQWakita,na.rm=TRUE),arr.ind=TRUE)[1,]
    }
    else if(q=="wakita2"){
      DeltaQWakita <- DeltaQ
      for(i in 1:length(DeltaQ[,1])){
        index <- which.max(DeltaQ[i,])
        DeltaQWakita[i,] <- NA
        DeltaQWakita[i,index]<-DeltaQ[i,index]
      }
      for(i in 1:(length(DeltaQ[,1])-1)){
        sizeCommI <- length(na.omit(DeltaQ[i,]))
        for(j in (i+1):length(DeltaQ[,1])){
          if(!is.na(DeltaQ[i,j])){
            sizeCommJ <- length(na.omit(DeltaQ[j,]))
            consol <- min(sizeCommI/sizeCommJ,sizeCommJ/sizeCommI)
            DeltaQWakita[i,j] <- DeltaQ[i,j]*consol
          }
        }  
      }
      maxEl <- which(DeltaQWakita==max(DeltaQWakita,na.rm=TRUE),arr.ind=TRUE)[1,]
    }
    else if(q=="wakita3"){
      DeltaQWakita <- DeltaQ
      for(i in 1:(length(DeltaQ[,1])-1)){
        sizeCommI <- length(which(C==i))
        for(j in (i+1):length(DeltaQ[,1])){
          if(!is.na(DeltaQ[i,j])){
            sizeCommJ <- length(which(C==j))
            consol <- min(sizeCommI/sizeCommJ,sizeCommJ/sizeCommI)
            DeltaQWakita[i,j] <- DeltaQ[i,j]*consol
          } 
        }  
      }
      maxEl <- which(DeltaQWakita==max(DeltaQWakita,na.rm=TRUE),arr.ind=TRUE)[1,]
    }
    
    if(randomized !=0){
      if(randomized>length(DeltaQ[,1])){
        randomized <- length(DeltaQ[,1])
      }
      rows <- sample(1:length(DeltaQ[,1]),randomized,replace=FALSE)
      DeltaQRandom <- DeltaQ[rows,]
      if(randomized>1){
        tmpEl <- which(DeltaQRandom==max(DeltaQRandom,na.rm=TRUE),arr.ind=TRUE)[1,]
        maxEl <- c(rows[tmpEl[1]],tmpEl[2])
      }
      else{
        tmpEl <- which(DeltaQRandom==max(DeltaQRandom,na.rm=TRUE))
        maxEl <- c(rows,tmpEl)
      }
    }
    
    maxI <- maxEl[1]
    maxJ <- maxEl[2]
    maxQ <- DeltaQ[maxI,maxJ]
    if(maxI>maxJ){
      tmp <- maxI
      maxI <- maxJ
      maxJ <- tmp
    }
    
    DeltaQ <- updateDeltaQ(DeltaQ,a,maxI,maxJ)     
    
    a[maxI] <- a[maxI]+a[maxJ]
    a <- a[-maxJ]
    
    C_join <- C
    indices <- which(C==maxJ)
    C_join[indices] <- maxI
    C_join <- updateC(C_join)
    
    Q_join <- Q+maxQ
    
    if(coarse!=0){
      newNumComm <- max(C_join)
      change <- numComm-newNumComm
      percentage <- change/numComm
      if(percentage>coarse){
        
        #build new network
        currA <- get.adjacency(network,attr="weight")
        adjacency <- matrix(0,ncol=max(C_join),nrow=max(C_join))
        for(i in 1:max(C_join)){
          verticesI <- which(C_join==i)
          for(j in 1:max(C_join)){
            verticesJ <- which(C_join==j)
            edge <- 0
            for(v in verticesI){
              for(u in verticesJ){
                if(are.connected(network,v,u)){
                  edge <- edge + currA[v,u]
                }
              }
            }
            adjacency[i,j] <- edge
            adjacency[j,i] <- adjacency[i,j]
          }
          adjacency[i,i] <- adjacency[i,i]/2
        }
        newGraph <- graph.adjacency(adjacency, mode="undirected", weighted=TRUE)
        tmpC <- initialSingletons(vcount(newGraph))
        tmpQ <- calculateQ(newGraph,tmpC)
        
        if(refine=="complete"){
          tmp <- completeGreedyRefinement(newGraph,tmpC,tmpQ)     
        }
        else if(refine=="fast"){
          tmp <- fastGreedyRefinement(newGraph,tmpC,tmpQ)
        }
        else if(refine=="kernighan"){
          tmp <- kernighanLinRefinement(newGraph,tmpC,tmpQ)
        }
        tmpC <- tmp[1:vcount(newGraph)]
        tmpQ <- tmp[vcount(newGraph)+1]
        
        #update the community structure for the original network
        C_tmp <- C_join
        for(c in 1:max(tmpC)){
          indices <- which(tmpC==c)
          vertices <- NULL
          for(i in indices){
            vertices <- c(vertices,which(C_join==i))
          }
          C_tmp[vertices] <- c
        }
        C_join <- updateC(C_tmp)
        #print(C_join)
        Q_join <- calculateQ(network,C_join)
        DeltaQ <- calculateDeltaQ(network,C_join)
        a <- calculateA(network,C_join)
        
        numComm <- newNumComm
      }
    }
    
    if(refine=="adclust"){
      tmp <- fastGreedyRefinement(network,C_join,Q_join)
      C_join <- tmp[1:vcount(network)]
      Q_join <- tmp[vcount(network)+1]
      DeltaQ <- calculateDeltaQ(network,C_join)
      a <- calculateA(network,C_join)
      
    }
    
    if(Q_join>Q_best){
      Q_best <- Q_join
      C_best <- C_join
    }
    
    Q <- Q_join
    C <- C_join     
  }
  
  if(refine=="complete"){
    tmp <- completeGreedyRefinement(network,C_best,Q_best)
    C_best <- tmp[1:vcount(network)]
    Q_best <- tmp[vcount(network)+1]
  }
  else if(refine=="fast"){
    tmp <- fastGreedyRefinement(network,C_best,Q_best)
    C_best <- tmp[1:vcount(network)]
    Q_best <- tmp[vcount(network)+1] 
  }
  else if(refine=="kernighan"){
    tmp <- kernighanLinRefinement(network,C_best,Q_best)
    C_best <- tmp[1:vcount(network)]
    Q_best <- tmp[vcount(network)+1]
  }
  
  result <- c(C_best,Q_best)
  return(result)
}

#RG+ approach using the randomized greedy approach
callRgplus <- function(network,z,randomized){
  
  #print(1)
  res <- callGreedy(network,randomized=randomized)
  C_cores <- res[1:vcount(network)]
  
  for(iter in 2:z){
    #print(iter)
    tmp <- callGreedy(network,randomized=randomized)
    C_tmp <- tmp[1:vcount(network)]
    C_new <- seq(0,0,length.out=vcount(network))
    for(i in 1:(vcount(network)-1)){
      for(j in (i+1):vcount(network)){
        if(C_cores[i]==C_cores[j]&C_tmp[i]==C_tmp[j]){
          C_new[i] <- C_cores[i]
          C_new[j] <- C_cores[i]
        }
      }  
    }
    zeros <- which(C_new==0)
    for(index in zeros){
      C_new[index] <- max(C_new)+1
    }
    
    C_cores <- updateC(C_new)
    
  }
  #print("final")
  result <- callGreedy(network,initialC=C_cores,randomized=randomized)
  return(result)
}

#Complete greedy refinement
completeGreedyRefinement <- function(network,C,Q){
  
  changed <- 1
  
  while(changed==1){
    changed <- 0
    DeltaQ_max <- -Inf
    i_max <- -1
    j_max <- -1
    for(i in 1:vcount(network)){
      neighs <- neighbors(network,i)
      connectedC <- C[neighs]
      connectedC <- connectedC[connectedC!=C[i]]
      connectedC <- unique(connectedC)
      for(j in connectedC){
        Qmove <- calculateQMove(network,C,j,i)
        if(Qmove>DeltaQ_max){
          DeltaQ_max <- Qmove
          i_max <- i
          j_max <- j
        }
      }
    }  
    
    if(DeltaQ_max>0){
      C[i_max]<-j_max
      C <- updateC(C)
      changed <- 1
      Q <- Q+DeltaQ_max
    }
  }
  result <- c(C,Q)
  return(result)
}

#Fast greedy refinement
fastGreedyRefinement <- function(network,C,Q){
  
  degrees <- degree(network,1:vcount(network))
  sortedDegrees <- sort(degrees,decreasing=FALSE)
  sortedDegrees <- unique(sortedDegrees)
  indices <- NULL
  
  for(s in sortedDegrees){
    indices <- c(indices, which(degrees==s))
  }
  
  changed<-1
  while(changed==1){
    changed <- 0
    for(i in indices){
      DeltaQ_max <- -Inf
      j_max <- -1
      neighs <- neighbors(network,i)
      connectedC <- C[neighs]
      connectedC <- connectedC[connectedC!=C[i]]
      connectedC <- unique(connectedC)
      for(j in connectedC){
        Qmove <- calculateQMove(network,C,j,i)
        if(Qmove>DeltaQ_max){
          DeltaQ_max <- Qmove
          j_max <- j
        }
      }
      if(DeltaQ_max>0){
        C[i] <- j_max
        C <- updateC(C)
        Q <- Q+DeltaQ_max
        changed <- 1
      }
    }
  }
  
  result <- c(C,Q)
  return(result)
}

#Kernighan Lin refinement
kernighanLinRefinement <- function(network,C,Q){
  
  changed <- 1
  while(changed==1){
    changed <- 0
    C_best <- C
    Q_best <- Q
    for(i in 1:vcount(network)){
      DeltaQ_max <- -Inf
      j_max <- -1
      neighs <- neighbors(network,i)
      connectedC <- C[neighs]
      connectedC <- connectedC[connectedC!=C[i]]
      connectedC <- unique(connectedC)
      for(j in connectedC){
        Qmove <- calculateQMove(network,C,j,i)
        if(Qmove>DeltaQ_max){
          DeltaQ_max <- Qmove
          j_max <- j
        }
      }
      C[i] <- j_max
      C <- updateC(C)
      Q <- Q + DeltaQ_max
      if(Q>Q_best&all.equal(Q,Q_best)==F){
        Q_best <- Q
        C_best <- C
        changed <- 1
      }
    }
    C <- C_best
    Q <- Q_best
  }
  
  result <- c(C,Q)
  return(result)
}

#MSG-VM algorithm
callMsgvm <- function(network,initialC=seq(0,0,length.out=vcount(network)), parL){
  
  if(initialC[1]==0){
    C <- initialSingletons(vcount(network))
    DeltaQ <- matrix(,nrow=max(C),ncol=max(C))
    a <- seq(0,0,length.out=max(C))
    
    for(i in 1:length(a)){
      a[i] <- getDegree(network,i)/(sum(get.adjacency(network,attr="weight")))  
    }
    
    for(i in 1:(max(C)-1)){
      for(j in (i+1):max(C)){
        if(are.connected(network,i,j)){
          DeltaQ[i,j] <- 2/(sum(get.adjacency(network,attr="weight")))-2*getDegree(network,i)*getDegree(network,j)/(sum(get.adjacency(network,attr="weight"))*sum(get.adjacency(network,attr="weight")))
          DeltaQ[j,i] <- DeltaQ[i,j]
        }
      }  
    }
  }
  else{
    C <- initialC
    
    a <- calculateA(network,C)
    DeltaQ <- calculateDeltaQ(network,C)
  }
  
  #print(DeltaQ)
  levelSet <- getLevelSet(DeltaQ)
  
  Q <- calculateQ(network,C)
  Q_best <- Q
  C_best <- C
  elements <- TRUE
  
  while(elements){
    if(length(levelSet)>3){
      if(levelSet[3,1]>0){
        #percentage <- (vcount(network)-max(C))/vcount(network)*100
        #print(paste(percentage,"%",sep=""))
        
        touched <- seq(0,0,length.out=max(C))
        if(length(levelSet[1,])<parL){
          parL <- length(levelSet[1,])
        }
        mp <- levelSet[,1:parL]
        mp <- mp[,mp[3,]>0]
        #print(levelSet)
        #print(mp)
        
        C_join <- C
        if(length(mp)>3){
          lengthMp <- length(mp[1,])
        }
        else{
          lengthMp <- 1
        }
        
        for(e in 1:lengthMp){
          if(lengthMp>1){
            i <- mp[1,e]
            j <- mp[2,e]
          }
          else{
            i <- mp[1]
            j <- mp[2]
          }  
          
          #print(length(touched))
          #print(i)
          #print(j)
          if(touched[i]==0&touched[j]==0){
            
            indexI <- min(which(C==i))
            i_join <- C_join[indexI]
            
            indexJ <- min(which(C==j))
            j_join <- C_join[indexJ]
            
            DeltaQ <- updateDeltaQ(DeltaQ,a,i_join,j_join)
            for(k in 1:length(DeltaQ[i_join,])){
              levelSet <- updateLevelSet(levelSet,i_join,k,DeltaQ[i_join,k])
            }
            
            #print(j_join)
            #print(levelSet)
            if(length(levelSet)>3){
              levelSet <- levelSet[,levelSet[1,]!=j_join]
              levelSet[1,] <- sapply(levelSet[1,], function(x) if(x>j_join){x=x-1}else{x=x})
              levelSet <- levelSet[,levelSet[2,]!=j_join]
              #print(levelSet)
              if(length(levelSet)>3){
                levelSet[2,] <- sapply(levelSet[2,], function(x) if(x>j_join){x=x-1}else{x=x})
              }
              else{
                if(levelSet[2]==j_join){
                  elements <- FALSE
                }
                else{
                  if(levelSet[2]>j_join){
                    levelSet[2] <- levelSet[2]-1
                  }
                  
                }
              }
            }  
            else{
              if(levelSet[1]==j_join|levelSet[2]==j_join){
                elements <- FALSE
              }
              else{
                if(levelSet[1]>j_join){
                  levelSet[1] <- levelSet[1]-1
                }
                if(levelSet[2]>j_join){
                  levelSet[2] <- levelSet[2]-1
                }
              }
            }
            
            a[i_join] <- a[i_join]+a[j_join]
            a <- a[-j_join]
            
            touched[i] <- 1
            touched[j] <- 1
            
            indices <- which(C==j)
            C_join[indices] <- i_join
            C_join <- updateC(C_join)  
          }
          
        }
      }
      else{
        break
      }
    }
    else{
      if(levelSet[3]>0){
        elements <- FALSE
        C_join <- sapply(C_join, function(x) x=1)
      }
      else{
        break
      }
    }

    C <- C_join
    Q <- calculateQ(network,C)
    
    if(Q>Q_best){
      Q_best <- Q
      C_best <- C
    }
  }
  
  result <- fastGreedyRefinement(network,C_best,Q_best)
  return(result)
}

#calculate the level set used in the msg-vm algorithm for a given matrix
getLevelSet <- function(matrix){
  levelSet <- NULL
  
  for(i in 1:(length(matrix[,1])-1)){
    #percentage <- (i-1)/length(matrix[,1]-1)*100
    #print(paste(percentage,"%",sep=""))
    for(j in (i+1):length(matrix[1,])){
      if(!is.na(matrix[i,j])){
        nextEl <- rbind(i,j,matrix[i,j])
      }
      else{
        nextEl <- rbind(i,j,0)
      }
      levelSet <- cbind(levelSet,nextEl)
    }
    
  }
  levelSet <- levelSet[,order(levelSet[3,],decreasing=T)]
  return(levelSet)
}

updateLevelSet <- function(levelSet,i,j,DeltaQ){
  
  if(i>j){
    tmp <- i
    i <- j
    j <- tmp
  }
  
  index <- which(levelSet[1,]==i&levelSet[2,]==j)
  if(!is.na(DeltaQ)){
    levelSet[3,index] <- DeltaQ
  }
  else{
    levelSet[3,index] <- 0
  }
  
  
  levelSet <- levelSet[,order(levelSet[3,],decreasing=T)]
  return(levelSet)
}

#CD algorithm
callCD <- function(network,initialC=seq(0,0,length.out=vcount(network)),maxC,iter,p){
  
  if(initialC[1]==0){
    if(maxC<vcount(network)){
      C <- sample(1:maxC,vcount(network),replace=TRUE)
    }
    else{
      C <- initialSingletons(vcount(network))
    }
  }
  else{
    C <- initialC
  }
  
  C_best <- C
  Q_best <- calculateQ(network,C)
  
  C_curr <- C_best
  Q_curr <- Q_best
  
  for(c in 1:iter){
    #percentage <- c/iter*100
    #print(paste(percentage,"%"))
    
    tmp <- completeGreedyRefinement(network,C_curr,Q_curr)
    C_curr <- tmp[1:vcount(network)]
    Q_curr <- tmp[vcount(network)+1]
    
    if(Q_curr>Q_best){
      C_best <- C_curr
      Q_best <- Q_curr
    }
    for(v in 1:vcount(network)){
      move <- sample(c(0,1),1,prob=c((1-p),p))
      if(move==1){
        community <- sample(1:(max(C_curr)+1),1)
        C_curr[v] <- community
        C_curr <- updateC(C_curr)
      } 
    }
    Q_curr <- calculateQ(network,C_curr)
  }
  
  result <- c(C_best,Q_best)
  return(result)
  
}

#greedy approach (louvain) according to Blondel et al.
callLouvain <- function(network,initialC=seq(0,0,length.out=vcount(network))){
  
  if(initialC[1]==0){
    C <- initialSingletons(vcount(network))
    
  }
  else{
    C <- initialC
  }
  
  Q <- calculateQ(network,C)
  currGraph <- network
  C_new <- C
  Q_new <- Q
  steps <- 1
  
  changed <- 1
  while(changed==1){
    changed <- 0
    tmp <- fastGreedyRefinement(currGraph,C_new,Q_new)
    C_curr <- tmp[1:vcount(currGraph)]
    Q_curr <- tmp[vcount(currGraph)+1]
    
    same <- 1
    for(i in 1:max(C_new)){
      vertInComm <- sort(which(C_new==i))
      newComm <- C_curr[vertInComm[1]]
      vertInNewComm <- sort(which(C_curr==newComm))
      bool <- vertInComm==vertInNewComm
      if(length(which(bool==FALSE))>0){
        same <- 0
        break
      }
    }
    
    #update the community structure for the original network
    C_tmp <- C
    for(c in 1:max(C_curr)){
      indices <- which(C_curr==c)
      vertices <- NULL
      for(i in indices){
        vertices <- c(vertices,which(C==i))
      }
      C_tmp[vertices] <- c
    }
    C <- updateC(C_tmp)
    
    if(same==0){
      #build new network
      currA <- get.adjacency(currGraph,attr="weight")
      adjacency <- matrix(0,ncol=max(C_curr),nrow=max(C_curr))
      for(i in 1:max(C_curr)){
        verticesI <- which(C_curr==i)
        for(j in 1:max(C_curr)){
          verticesJ <- which(C_curr==j)
          edge <- 0
          for(v in verticesI){
            for(u in verticesJ){
              if(are.connected(currGraph,v,u)){
                edge <- edge + currA[v,u]
              }
            }
          }
          adjacency[i,j] <- edge
          adjacency[j,i] <- adjacency[i,j]
        }
        adjacency[i,i] <- adjacency[i,i]/2
      }
      newGraph <- graph.adjacency(adjacency, mode="undirected", weighted=TRUE)
      currGraph <- newGraph
      
      C_new <- initialSingletons(vcount(newGraph))
      Q_new <- calculateQ(newGraph,C_new)
      
      changed <- 1
    }
  }
  
  Q <- calculateQ(network,C)
  result <- c(C,Q)
  return(result)
  
}

#greedy appraoch based on vertex similiarty and using a hybrid merge
callVertexSim <- function(network,frac){
  weightedNetwork <- weighting(network)
  
  #identification of perliminary communities
  C <- seq(0,0,length.out=vcount(network))
  commNum <- 1
  weightes <- weightedNetwork[,order(weightedNetwork[3,],decreasing=TRUE)]
  
  for(w in 1:length(weightes[1,])){
    v <- weightes[1,w]
    u <- weightes[2,w]
    
    if(C[v]==0& C[u]==0){
      C[v] <- commNum
      C[u] <- commNum
      commNum <- commNum+1
    }
  }
  missing <- which(C==0)
  for(m in missing){
    C[m] <- commNum
    commNum <- commNum+1
  }
  
  #merge communities
  Q <- calculateQ(network,C)
  C_best <- C
  Q_best <- Q
  
  a <- calculateA(network,C)
  DeltaQ <- calculateDeltaQ(network,C)
  
  for(iteration in 1:ceiling(log(vcount(network)))){
    #percentage <- (iteration-1)/ceiling(log(vcount(network)))*100
    #print(paste(percentage,"%",sep=""))
    
    C_join <- C
    linkSet <- matrix(,nrow=max(C),ncol=max(C))
    
    for(i in 1:max(C)){
      #percent <- (i-1)/max(C)*100
      #print(paste(percent,"%",sep=""))
      #print(max(C))
      if(max(C)>1){
        j_max <- which(DeltaQ[i,]==max(DeltaQ[i,],na.rm=TRUE))[1]
        DeltaQ_max <- DeltaQ[i,j_max]
        
        if(DeltaQ_max>0){
          linkSet[i,j_max] <- 1
        }
      }
    }
    
    pairwise <- frac*ceiling(log(vcount(network)))
    if(iteration<=pairwise){
      for(i in 1:(max(C)-1)){
        for(j in (i+1):max(C)){
          if(!is.na(linkSet[i,j])&!is.na(linkSet[j,i])){
            
            indexI <- min(which(C==i))
            i_join <- C_join[indexI]
            
            indexJ <- min(which(C==j))
            j_join <- C_join[indexJ]
            
            verticesJ <- which(C==j)
            C_join[verticesJ] <- i_join
            C_join <- updateC(C_join)
            Q <- calculateQ(network,C_join)
            DeltaQ <- updateDeltaQ(DeltaQ,a,i_join,j_join)
            a[i_join] <- a[i_join]+a[j_join]
            a <- a[-j_join]
            
            if(Q>Q_best){
              Q_best <- Q
              C_best <- C_join
            }
            linkSet[i,j] <- NA
            linkSet[j,i] <- NA
          }
        }
      }
      
    }else{
      for(i in 1:max(C)){
        for(j in 1:max(C)){
          if(!is.na(linkSet[i,j])){
            if(i>j){
              tmp <- i
              i <- j
              j <- tmp
            }
            indexI <- min(which(C==i))
            i_join <- C_join[indexI]
            
            indexJ <- min(which(C==j))
            j_join <- C_join[indexJ]
            verticesJ <- which(C_join==j_join)
            C_join[verticesJ] <- i_join
            C_join <- updateC(C_join)
            Q <- calculateQ(network,C_join)
            DeltaQ <- updateDeltaQ(DeltaQ,a,i_join,j_join)
            a[i_join] <- a[i_join]+a[j_join]
            a <- a[-j_join]
            
            if(Q>Q_best){
              Q_best <- Q
              
              C_best <- C_join
            }
            linkSet[i,j] <- NA
            linkSet[j,i] <- NA
          }
        }
      }
    }
    C <- C_join   
  }
  
  result <- c(C_best,Q_best)
  return(result)
  
}

weighting <- function(network){
  
  edges <- get.edgelist(network)
  weightedEdges <- NULL
  for(e in 1:length(edges[,1])){
    v1 <- edges[e,1]
    #v1 <- unlist(strsplit(v1,split=""))[2]
    #v1 <- as.numeric(v1)
    v2 <- edges[e,2]
    #v2 <- unlist(strsplit(v2,split=""))[2]
    #v2 <- as.numeric(v2)
    col <- rbind(v1,v2,0)
    weightedEdges <- cbind(weightedEdges,col)
  }
  
  for(k in 1:ceiling(log(vcount(network)))){
    tab <- matrix(0,nrow=2,ncol=vcount(network))
    for(i in 1:vcount(network)){
      neighs <- neighbors(network,i)
      for(j in neighs){
        tab[2,j] <- i
      }
      for(j in neighs){
        col <- which((weightedEdges[1,]==i&weightedEdges[2,]==j)|(weightedEdges[1,]==j & weightedEdges[2,]==i))
        
        if(isTRUE(tab[1,j]==0 & weightedEdges[3,col]==0)){
          cFriends <- 0
          neighsJ <- neighbors(network,j)
          for(z in neighsJ){
            if(tab[2,z]==i){
              cFriends <- cFriends+1
            }
          }
          cosine <- cFriends/(sqrt(length(neighs)*length(neighsJ)))
          weightedEdges[3,col] <- cosine
          tab[1,j]<-1
        }
      }
    }
  }
  return(weightedEdges)
}

#mome approach
callMome <- function(network){
  
  adjacency <- get.adjacency(network,attr="weight")
  communities <- NULL
  res <- computeCoarsening(adjacency,adjacency,communities)
  coarsening <- res[1:(length(res)/2)]
  communities <- res[(length(res)/2+1):length(res)]
  
  C_tmp <- 0
  
  for(g in length(coarsening):2){
    #percentage <- (length(coarsening)-g)/length(coarsening)*100
    #print(paste(percentage,"%",sep=""))
    
    graph <- unlist(coarsening[g])
    dim <- sqrt(length(graph))
    adj <- matrix(graph,nrow=dim,byrow=TRUE)
    currGraph <- graph.adjacency(adj,mode="undirected",weighted=TRUE)
    C_curr <- unlist(communities[g])
    if(C_tmp[1] !=0){
      for(c in 1:max(C_tmp)){
        indices <- which(C_tmp==c)
        vertices <- NULL
        for(i in indices){
          vertices <- c(vertices,which(C_curr==i))
        }
        C_curr[vertices] <- c
      }
    }
    C_curr <- updateC(C_curr)
    Q_curr <- calculateQ(currGraph,C_curr)
    tmp <- fastGreedyRefinement(currGraph,C_curr,Q_curr)
    C_tmp <- tmp[1:vcount(currGraph)]
  }
  
  C <- unlist(communities[1])
  for(c in 1:max(C_tmp)){
    indices <- which(C_tmp==c)
    vertices <- NULL
    for(i in indices){
      vertices <- c(vertices,which(C==i))
    }
    C[vertices] <- c
  }
  Q <- calculateQ(network,C)
  
  result <- c(C,Q)
  return(result)
}

computeCoarsening <- function(currA,coarsening,communities){
  
  #print("Coarsening...")
  currGraph <- graph.adjacency(currA, mode="undirected", weighted=TRUE)
  C <- initialSingletons(vcount(currGraph))
  
  for(i in 1:vcount(currGraph)){
    Q_max <- -Inf
    j_max <- -1
    commI <- C[i]
    neighs <- neighbors(currGraph,i)
    for(j in neighs){
      commJ <- C[j]
      if(commJ!=commI){
        Q_merge <- calculateQMerge(currGraph,C,commI,commJ)
        if(Q_merge>Q_max){
          Q_max <- Q_merge
          j_max <- j
        }
      }
    }
    if(Q_max>0){
      commJ <- C[j_max]
      verticesJ <- which(C==commJ)
      C[verticesJ] <- C[i]  
      C <- updateC(C)
    }    
  }
  
  #build new network
  adjacency <- matrix(0,ncol=max(C),nrow=max(C))
  for(i in 1:max(C)){
    verticesI <- which(C==i)
    for(j in 1:max(C)){
      verticesJ <- which(C==j)
      edge <- 0
      for(v in verticesI){
        for(u in verticesJ){
          if(are.connected(currGraph,v,u)){
            edge <- edge + currA[v,u]
          }
        }
      }
      adjacency[i,j] <- edge
      adjacency[j,i] <- adjacency[i,j]
    }
    adjacency[i,i] <- adjacency[i,i]/2
  }
  communities <- append(communities,list(C))
  
  if(length(adjacency[1,])<vcount(currGraph)){
    coarsening <- append(coarsening,list(adjacency))
    result <- computeCoarsening(adjacency,coarsening,communities)
  }
  else{
    result <- c(coarsening,communities)
  }
  
  return(result)
}
