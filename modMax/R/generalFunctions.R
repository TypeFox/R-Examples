##definition of different initial partitions##
#place all vertices in their own communities
initialSingletons <- function(n){
  C <- NULL
  for(i in 1:n){
    C[i] <- i
  }
  return(C)
}

#initialze using prior knowledge according to Du et al.
initialPriorKnowledge <- function(network){
  C <- initialSingletons(vcount(network))
  currGraph <- network
  vertices <- c(1:vcount(network))
  edges <- get.edgelist(network)
  
  while(length(edges)>1){
    #percentage <- (vcount(network)-max(C))/vcount(network)*100
    #print(paste(percentage,"%",sep=""))
    
    degrees <- degree(currGraph,vertices)
    index <- which.max(degrees)
    vertex <- vertices[index]

    comm <- C[vertex]
    neighs <- neighbors(currGraph,vertex)
    C[neighs] <- comm
    C <- updateC(C)
    vertices <- vertices[vertices != vertex]
    if(length(edges)>2){
      deleteEdges <- which(edges==vertex,arr.ind=TRUE)
      if(length(deleteEdges)>2){
        deleteEdges <- deleteEdges[,1]
      }
      else{
        deleteEdges <- deleteEdges[1]
      }
      edges <- edges[-deleteEdges,]
    }
    else{
      if(edges[1]==vertex|edges[2]==vertex){
        edges <- NULL
        deleteEdges <- 1
      }      
    }
    currGraph <- delete.edges(currGraph,deleteEdges)
    for(n in neighs){
      vertices <- vertices[vertices != n]
      #print(edges)
      #print(n)
      if(degree(currGraph,n)>0){
        if(length(edges)>2){
          deleteEdges <- which(edges==n,arr.ind=TRUE)
          if(length(deleteEdges)>2){
            deleteEdges <- deleteEdges[,1]
          }
          else{
            deleteEdges <- deleteEdges[1]
          }
          edges <- edges[-deleteEdges,]
        }
        else{
          if(edges[1]==n|edges[2]==n){
            edges <- NULL
            deleteEdges <- 1
          } 
        }
        #print(deleteEdges)
        currGraph <- delete.edges(currGraph,deleteEdges)
      }
    }
  }
  return(C)
}

#initialze using random walkers according to Pujol et al.
initialRandomWalkers <- function(network){
  parT <- 3
  R <- 0.2
  
  adjacency <- get.adjacency(network,attr="weight")
  identity <- diag(vcount(network))
  D <- matrix(0,nrow=vcount(network),ncol=vcount(network))
  for(i in 1:vcount(network)){
    D[i,i] <- 1+ getDegree(network,i)
  }
  M <- as.matrix((adjacency+identity)%*%solve(D))
  degrees <- degree.distribution(network)[2:length(degree.distribution(network))]
  index <- length(degrees)
  usedDeg <- NULL
  thresh <- 0
  while(thresh+degrees[index]<=R){
    thresh <- thresh+degrees[index]
    index <- index-1
  }

  seeds <- NULL
  for(i in 1:vcount(network)){
      if(degree(network,i)>=index){
        seeds <- c(seeds,i)
      }
  }
  G_old <- matrix(0,nrow=vcount(network),ncol=length(seeds))
  j <- 1
  for(s in seeds){
    G_old[s,j] <- 1
    j <- j+1
  }

  for(t in 1:parT){
   G_new <- t(M)%*%G_old
   G_old <- G_new
  }
 C <- seq(0,0,length.out=vcount(network))
  for(i in 1:vcount(network)){
   comm <- which.max(G_new[i,])
   C[i] <- comm
  }
  commNum <- max(C)+1
  vertices <- which(C==0)
  for(v in vertices){
   C[v] <- commNum
   CommNum <- commNum+1
  }
  C <- updateC(C)
  
  return(C)  
}

#initialize using subgraph similarity according to Xiang et al.
initialSubgraphSim <- function(network){
  
  C <- initialSingletons(vcount(network))
  merges <- seq(0,0,length.out=vcount(network))
  adjacency <- get.adjacency(network)
  
  for(i in 1:max(C)){
    #percentage <- (i-1)/max(C)*100
    #print(paste(percentage,"%",sep=""))
    s_max <- 0
    j_max <- -1
    neighs <- neighbors(network,i)
    for(j in 1:max(C)){
      if(j!=i){
        neighsJ <- neighbors(network,j)
        commonNeigh <- 0
        
        for(n in neighs){
          if(n %in% neighsJ){
            commonNeigh <- commonNeigh+1  
          }
        }
        s <- (adjacency[i,j]+commonNeigh)/sqrt(degree(network,i)*degree(network,j))
        if(s>s_max){
          s_max <- s
          j_max <- j
        }
      }
    }
    merges[i] <- j_max
  }

  commNum <- 1
  
  for(i in 1:vcount(network)){
    
    if(i<merges[i]){
      C[merges[i]] <- C[i]
    }
    else{
      C[i] <- C[merges[i]]
    }
  }
  C <- updateC(C)
  
  return(C)
}

#randomly determine an initial partition
initialRandom <- function(network){
  C <- 1
  while(max(C)==1 | min(C)==vcount(network)){
    C <- sample(1:vcount(network),vcount(network), replace=TRUE)
  }
  
  return(C)
}

#use community structure from the greedy algorithm as initial partition
initialGreedy <- function(network){
  result <- callGreedy(network)
  C <- result[1:vcount(network)]
  return(C)
}

#place all vertices in one community
initialWholeGraph <- function(network){
  C <- NULL
  for(i in 1:vcount(network)){
    C[i] <- 1
  }
  return(C)
}

#calculate the elements of vector a
calculateA <- function(network,C){
  
  a <- seq(0,0,length.out=max(C))
  for(i in 1:length(a)){
    vertices <- which(C==i)
    for(v in vertices){
      a[i] <- a[i]+getDegree(network,v)
    }
    a[i] <- a[i]/(sum(get.adjacency(network,attr="weight")))
  }
  
  return(a)
}

#calculate the matrix Delta Q
calculateDeltaQ <- function(network,C){
  #print("calculate Q...")
  DeltaQ <- matrix(,nrow=max(C),ncol=max(C))
  
  if(max(C)>1){
    for(i in 1:(max(C)-1)){
      
      for(j in (i+1):max(C)){
        value <- calculateQMerge(network,C,i,j)
        if(value!=0){
          DeltaQ[i,j] <- value
          DeltaQ[j,i] <- DeltaQ[i,j]
        }      
      }
      #percentage <- i/(max(C)-1)*100
      #print(paste(percentage,"%",sep=""))
    }
  }
  
  return(DeltaQ)
}

#general calculation of Q
calculateQ <- function(network,C){
  
  A <- get.adjacency(network,attr="weight")
  Q <- 0
  
  for(c in 1:max(C)){
    a <- 0
    e <- 0
    for(i in 1:vcount(network)){
      if(C[i]==c){
        a <- a+getDegree(network,i)
        neighbours <- neighbors(network,i)
        if(length(neighbours)>0){
          for(j in 1:length(neighbours)){
            if(C[neighbours[j]]==c){
              e <- e+A[i,as.numeric(neighbours[j])]
            }
          }
        }
      }
    }

    a <- a/(sum(A))
    e <- e/(sum(A))
    Q <- Q+(e-a^2)
    #Q <- round(Q,2)
  }
  
  return(Q)
}

#calculate the change of modularity due to the merge of 2 communities
calculateQMerge <- function(network,C,id1,id2){
  
  adj <- get.adjacency(network,attr="weight")
  
  comm1 <- NULL
  comm2 <- NULL
  
  Q_merge <- 0
  
  for(i in 1:vcount(network)){
    if(C[i]==id1){
      comm1 <- cbind(comm1,i)
    } else if(C[i]==id2){
      comm2 <- cbind(comm2,i)
    }
  }
  
  e <-0
  connect <- 0
  for(i in 1:length(comm1)){
    for(j in 1:length(comm2)){
      if(are.connected(network,comm1[i],comm2[j])){
        #print(comm1[i])
        #print(comm2[j])
        e <- e+adj[comm1[i],comm2[j]]
        connect <- 1
      }
    }
  }
  if(connect==1){
    a_1 <-0
    for(i in 1:length(comm1)){
      a_1 <- a_1+getDegree(network,comm1[i])
    }
    a_2 <- 0
    for(i in 1:length(comm2)){
      a_2 <- a_2+getDegree(network,comm2[i])
    }
    
    e <- e/(sum(adj))
    a_1 <- a_1/(sum(adj))
    a_2 <- a_2/(sum(adj))
    
    Q_merge <- 2*(e-a_1*a_2)
  }

  return(Q_merge)
}

#calculate the change of modularity due to a move of a vertex to another community
calculateQMove <- function(network,C,community,vertex){
  
  adjacency <- get.adjacency(network,attr="weight")

  Q_move <- -Inf
  
  if(C[vertex]!= community){
    newComm <- NULL
    oldComm <- NULL
    for(i in 1:length(C)){
      if(C[i]==community){
        newComm <- cbind(newComm, i)  
      }
      if(C[i]==C[vertex]){
        oldComm <- cbind(oldComm, i)
      }
    }
    
    links_new <- 0
    if(length(newComm)>0){
      for(j in 1:length(newComm)){
        if(are.connected(network,vertex,newComm[j])){
          links_new <- links_new+adjacency[vertex,newComm[j]]
        }
      }
    }
    
    if(community==(max(C)+1) | links_new !=0){
      links_old <- 0
      for(i in 1:length(oldComm)){
        if(are.connected(network,vertex,oldComm[i])){
          links_old <- links_old+adjacency[vertex,oldComm[i]]
        }
      }
      a_old <- 0
      for(i in 1:length(oldComm)){
        a_old <- a_old+getDegree(network,oldComm[i])
      }
      
      a_new <- 0
      if(length(newComm)>0){
        for(j in 1:length(newComm)){
          a_new <- a_new+getDegree(network,newComm[j])
        }
      }
      
      Q_move <- 2*(links_new-links_old)/sum(adjacency)-2*(getDegree(network,vertex)*(a_new-(a_old-getDegree(network,vertex))))/(sum(adjacency)*sum(adjacency))
      #print(Q_move)
    }
  }
  else{
    Q_move <- 0
  }
  
  #print(Q_move)
  return(Q_move)
}

#update DeltaQ according to generic greedy algorithm
updateDeltaQ <- function(DeltaQ,a,maxI,maxJ){
  
  for(k in 1:length(DeltaQ[maxI,])){
    if(!is.na(DeltaQ[maxI,k]) & !is.na(DeltaQ[maxJ,k])){
      DeltaQ[maxI,k] <- DeltaQ[maxI,k]+DeltaQ[maxJ,k]
    }
    else if(!is.na(DeltaQ[maxI,k]) & is.na(DeltaQ[maxJ,k])){
      DeltaQ[maxI,k] <- DeltaQ[maxI,k]-2*(a[maxJ]*a[k])
    }
    else if(is.na(DeltaQ[maxI,k]) & !is.na(DeltaQ[maxJ,k])){
      DeltaQ[maxI,k] <- DeltaQ[maxJ,k]-2*(a[maxI]*a[k])
    }

    DeltaQ[k,maxI] <- DeltaQ[maxI,k]
  }
  
  DeltaQ[maxI,maxI] <- NA
  DeltaQ <- DeltaQ[-maxJ,]

  if(length(DeltaQ)>2){
    DeltaQ <- DeltaQ[,-maxJ]
  }
  else{
    DeltaQ <- DeltaQ[-maxJ]
  }
  return(DeltaQ)
}

#update the community ids so that they always start with 1,...
updateC <- function(C){
  
  currComm <- 1
  while(currComm < max(C)){
    if(currComm %in% C){
      currComm <- currComm+1
    } else{
      C <- sapply(C, function(x) if(x>currComm){x=x-1}else{x=x})
    }
  }
  return(C)
}

#calculate the modularity matrix B
calculateModMatrix<-function(network){
  A <- get.adjacency(network,attr="weight")
  B <- matrix(,nrow=length(A[,1]),ncol=length(A[1,]))
  for(i in 1:length(B[,1])){
    for(j in i:length(B[,1])){
      B[i,j] <- A[i,j]-getDegree(network,i)*getDegree(network,j)/(sum(A))
      B[j,i] <- B[i,j]
    }
    #percentage <- i/(length(B[,1]))*100
    #print(paste(percentage,"%",sep=""))
  }
  
  return(B)
}

#calculate random graph with same degree distribution as given graph
calculateRandomGraph <- function(network){
  
  degDis <- degree.distribution(network)
  degDis <- sapply(degDis,function(x) x=x*vcount(network))
  degDis <- degDis[-1]
  weight <- get.edge.attribute(network, "weight")
  
  degs <- NULL
  
  for(i in 1:length(degDis)){
    degs <- c(degs,seq(i,i,length.out=degDis[i]))
  }

  vertDegs <- sample(degs,vcount(network),replace=FALSE)
  rndNetwork <- degree.sequence.game(out.deg=vertDegs, method="vl")
  rndWeight <- sample(weight,ecount(rndNetwork),replace=FALSE)
  rndNetwork <- set.edge.attribute(rndNetwork,"weight",value=rndWeight)
  
  return(rndNetwork)
}

getDegree <- function(network,i){
  adj <- get.adjacency(network,attr="weight")
  deg <- sum(adj[i,])
  return(deg)
}

generateOutput <- function(res,rndRes=NULL,random=FALSE){
  
  result <- list()
  result <- append(result,list(max(res)))
  result <- append(result,list(res[length(res)]))
  if(random==TRUE){
    avg_random <- mean(rndRes)
    sd_random <- sd(rndRes)
    result <- append(result,list(avg_random))
    result <- append(result,list(sd_random))
  }
  result <- append(result,list(res[1:(length(res)-1)]))
  if(random==TRUE){
    result <- append(result,list(rndRes))
    names(result) <- c("number of communities","modularity","mean", 
                       "standard deviation", "community structure", 
                       "random modularity values")
  }
  else{
    names(result) <- c("number of communities","modularity", "community structure")
  }
  
  return(result)
}

