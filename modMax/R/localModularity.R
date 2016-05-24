#exploring of community structure around a certain vertex
localModularity <- function(adjacency,srcV,k){
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  res <- callLocalModularity(network,srcV,k)
  res[length(res)] <- round(res[length(res)],2)
  
  result <- generateLocalOutput(res)
  
  return(result)
}

#usage of local modularity to explore whole structure according to Wang et al.
localModularityWang <- function(adjacency,numRandom=0){
  
  network <- graph.adjacency(adjacency, mode="undirected",weighted=TRUE)
  res <- callLocalModularityWang(network)
  res[vcount(network)+1] <- round(as.numeric(res[vcount(network)+1]),2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callLocalModularityWang(rndNetwork)
      rndResult[vcount(rndNetwork)+1] <- round(as.numeric(rndResult[vcount(rndNetwork)+1]),2)
      randomResults <- rbind(randomResults,rndResult[vcount(rndNetwork)+1])
    }
    result <- generateOutput(res,randomResults,random=TRUE)
  }
  else{
    result <- generateOutput(res)
  }
  
  return(result)
  
}

callLocalModularity <- function(network,srcV,k){
  
  adjacency <- get.adjacency(network,attr="weight")
  C_l <- srcV
  B <- srcV
  neighbors <- neighbors(network,srcV)
  U <- neighbors
  R <- 0
  
  while(length(C_l)<k&(length(C_l)+length(U)<vcount(network))){
    #calculate T
    T <- 0
    edgesInB <- 0
    for(b in B){
      T <- T + getDegree(network,b)
      for(b2 in B){
        if(are.connected(network,b,b2)){
          edgesInB <- edgesInB +adjacency[b,b2]
        }
      }
    }
    T <- (T - edgesInB/2)*2
    deltaR_max <- -Inf
    j_max <- -1
    for(j in U){
      if(j != srcV){
        #calculate x
        x <- 0
        for(b in B){
          if(are.connected(network,j,b)){
            x <- x+adjacency[j,b]
          }
        }
        neighbors_j <- neighbors(network,j)
        neighbors_j <- neighbors_j[neighbors_j!=j]
        #calculate y and z
        jInB <- 0
        y <- 0
        z <- 0
        for(n in neighbors_j){
          if(!n %in% B){
            jInB <- 1
          }
          else if(n %in% B){
            inB <- 0
            for(u in U){
              if(u != j & are.connected(network,u,n)){
                inB <- 1
              }
            }
            if(inB==0){
              z <- z+getDegree(network,n)
            }
          }
        }
        if(jInB==1){
          y <- getDegree(network,j)
        }

        deltaR <- (x-R*y-z*(1-R))/(T-z+y)    
        if(deltaR > deltaR_max){
          deltaR_max <- deltaR
          j_max <- j
        }
      }
    }
    #update C_l
    C_l <- cbind(C_l,j_max)
    #update U
    neighbors_j <- neighbors(network,j_max)
    neighbors_j <- neighbors_j[neighbors_j!=j_max]
    U <- c(U,neighbors_j)
    U <- U[!U %in% C_l]
    U <- unique(U)
    #update B
    B <- NULL
    for(c in C_l){
      for(u in U){
        if(are.connected(network,c,u)){
          B <- cbind(B,c)
          break
        }
      }
    }
    R <- R + deltaR_max  
  }
  
  result <- c(C_l,R)
  return(result)
}

callLocalModularityWang <- function(network){
  tab <- NULL
  for(i in 1:vcount(network)){
    neighs <- neighbors(network,i)
    neighs <- neighs[neighs!=i]
    neighs <- paste(neighs,collapse=",")
    row <- c(i,neighs,getDegree(network,i),0)
    tab <- rbind(tab,row)
  }

  commNum <- 1
  while(0 %in% tab[,4]){
    #percentage <- length(which(tab[,4]!=0))/vcount(network)*100
    #print(paste(percentage,"%",sep=""))
    
    vertices <- which(tab[,4]==0)
    sourceV <- sample(vertices,1)
    tab <- localCommunityDetection(network,sourceV,commNum,tab)
    commNum <- commNum+1
  }
  C <- NULL
  for(i in 1:vcount(network)){
    C <- cbind(C,tab[i,4])
  }
  
  Q <- calculateQ(network,C)
  
  result <- c(C,Q)
  
  return(result)
}

localCommunityDetection<-function(network,sourceV,commNum,tab){
  
  Q_l <- 0
  tab[sourceV,4]<-commNum
  tab <- cbind(tab,0,0)
  S <- NULL
  neighs <- tab[sourceV,2]
  neighs <- as.numeric(unlist(strsplit(neighs,split=",")))
  for(i in neighs){
    if(tab[i,4]==0){
      S <- c(S,i)
      tab[i,6]<-1   #vertex has been visited
    }
  }
  
  adjacency <- get.adjacency(network,attr="weight")
  changed <- 1
  while(changed==1){
    unvisited <- which(tab[,6]==0)
    if(0 %in% tab[unvisited,4]){
      Q_max <- -Inf
      j_max <- -1
      for(j in S){
        #calculate L_in
        l_in <- 0
        indices <- which(tab[,4]==commNum)
        indices <- cbind(indices,j)
        for(v in indices){
          for(u in indices){
            if(are.connected(network,u,v)){
              l_in <- l_in+adjacency[u,v]
            }
          }
        }
        l_in <- l_in/2
        #calculate denominator
        deNom <- sum(as.numeric(tab[indices,3]))-l_in
        Q_new <- l_in/deNom
        if(Q_new < Q_l){
          tab[j,5]=as.numeric(tab[j,5])+1
          if(as.numeric(tab[j,5])==2){
            S <- S[S!=j]
          }
        }
        else{
          tab[j,5]=0
          if(Q_new > Q_max){
            Q_max <- Q_new
            j_max <- j
          }
        }
      }
      if(Q_max>Q_l){
        Q_l <- Q_max
        tab[j_max,4]=commNum
        #update S
        neighs <- tab[j_max,2]
        neighs <- as.numeric(unlist(strsplit(neighs,split=",")))
        for(i in neighs){
          if((as.numeric(tab[i,4])==0) & (as.numeric(tab[i,5])!=2)){
            S <- c(S,i)
            tab[i,6]<-1   #vertex has been visited
          }
        }
        S <- unique(S)
      }
      else{
        changed <- 0
      }
    }
    else{
      for(s in S){
        tab[s,4]=commNum
        changed <- 0
      }
    }
  }
  
  tab <- tab[,1:4]
  return(tab)
}

#generate special output for a local community structure
generateLocalOutput <- function(res){
  result <- list()
  result <- append(result, res[1:length(res)-1])
  result <- append(result, res[length(res)])
  names(result) <- c("Local community structure", "Local modularity")
  
  return(result)
}