#performance of individual node movement
individualNodeMovement <- function(network,C,Q,beta){
  
  vertex <- sample(1:vcount(network),1)
  community <- sample(1:(max(C)+1),1)
  accept <- 0  
  
  if(community != C[vertex]){
    
    #print(C)
    #print(vertex)
    #print(community)
    
    Q_move <- calculateQMove(network,C,community,vertex)
    #print(Q_move)
    connected <- 1
    C_tmp <- C
    C_tmp[vertex] <- community
    C_tmp <- updateC(C_tmp)
    
    for(i in 1:vcount(network)){
      connected_i <- 0
      neighs <- neighbors(network,i)
      for(j in neighs){
        if(C_tmp[j]==C_tmp[i]){
          connected_i <- 1
        }
      }
      if(connected_i==0){
        connected <- 0
        break
      }
    }

    if(connected==0){
      accept <- 0
    }
    else if(Q_move != -Inf){
      if(Q_move > 0){
        accept <- 1
      }
      else{
        prob <- exp(beta*Q_move)
        accept <- sample(c(0,1),1,replace=TRUE,c((1-prob),prob))
      }
    }
  }
  
  if(accept==1){
    C[vertex]= community
    
    Q <- Q + Q_move
    C <- updateC(C)
  }  
  result <- c(C,Q,accept) 
}

####performance of collective movements
#merging of communities
mergeCommunities <- function(network,C,Q,beta){
  
  accept <- 0
  if(max(C)>1){
    ids <- sample(1:max(C),2,replace=FALSE)
    Q_merge <- calculateQMerge(network,C,ids[1],ids[2])
  }
  else{
    Q_merge <- 0
  }
  
  #print(Q_merge)
  #print(Q_merge > 0)
  #print(all.equal(Q_merge,0))
  
  connected <- 1
  C_tmp <- C
  for(i in 1:length(C_tmp)){
    if(C_tmp[i]==ids[2]){
      C_tmp[i]==ids[1]
    }
  }
  C_tmp <- updateC(C_tmp)
  
  for(i in 1:vcount(network)){
    connected_i <- 0
    neighs <- neighbors(network,i)
    for(j in neighs){
      if(C_tmp[j]==C_tmp[i]){
        connected_i <- 1
      }
    }
    if(connected_i==0){
      connected <- 0
      break
    }
  }
  
  if(connected==0){
    accept <- 0
  }
  if(all.equal(Q_merge,0)==T){
    accept <- 0
  }
  else if(Q_merge > 0){
    accept <- 1
  }
  else{
    prob <- exp(beta*Q_merge)
    accept <- sample(c(0,1),1,replace=TRUE,c((1-prob),prob))
  }
  
  
  if(accept==1){
    for(i in 1:length(C)){
      if(C[i]==ids[2]){
        C[i]==ids[1]
      }
    }
    Q <- Q+Q_merge
    C <- updateC(C)
  }
  
  result <- c(C,Q,accept)
  return(result)
}

#splitting of communities
splitCommunities <- function(network,C,Q,beta,alpha,beta_sys){
  id <- sample(1:max(C),1)
  comm <- NULL
  for(i in 1:length(C)){
    if(C[i]==id){
      comm <- cbind(comm, i)
    }
  }
  
  if(length(comm)>1){
    sub <- induced.subgraph(network,comm)
    #plot(sub)
    C_sub <- sample(c(1,2),vcount(sub), replace=TRUE)
    C_sub <- updateC(C_sub)
    Q_sub <- calculateQ(sub,C_sub)
    
    while(beta<beta_sys){
      #print(Q_sub)
      result <- individualNodeMovement(sub,C_sub,Q_sub,beta)
      Q_sub <- result[vcount(sub)+1]
      #print(Q_sub)
      C_sub <- result[1:vcount(sub)]
      beta <- beta*alpha
    }
    
    Q_split <- Q_sub - Q
    
    accept <-0
    connected <- 1
    
    for(i in 1:vcount(sub)){
        connected_i <- 0
        neighs <- neighbors(sub,i)
        for(j in neighs){
          if(C_sub[j]==C_sub[i]){
            connected_i <- 1
          }
        }
        if(connected_i==0){
          connected <- 0
          break
        }
    }
    
    if(connected==0){
      accept <- 0
    }
    else if(Q_split>0){
      accept <- 1
    } else{
      prob <- exp(beta*Q_split)
      accept <- sample(c(0,1),1,replace=TRUE,c((1-prob),prob))
    }
    
    if(accept==1){
      newComm <- max(C)+1
      Q <- Q_sub
      for(i in 1:length(C_sub)){
        vertex=comm[i]
        if(C_sub[i]==1){
          C[vertex]=id
        }
        else{
          C[vertex]<- newComm
        }
        
      }
    }
    C <- updateC(C)
  }
  else{
    accept <- 0
  }

  result <- c(C,Q,accept)
  return(result)
}