#execute the generic extremal optimization algorithm
extremalOptimization <- function(adjacency, numRandom=0,
                                 refine=c("none","agents"),
                                 tau=FALSE,alpha_max=length(adjacency[1,]),steps=3){
  refine <- match.arg(refine)

  network <- graph.adjacency(adjacency,mode="undirected", weighted=TRUE)
  res <- callExtremalOptimization(network,refine=refine,tau=tau,alpha_max=alpha_max,steps)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  randomResults <- NULL
  if(numRandom!=0){
    for(i in 1:numRandom){
      rndNetwork <- calculateRandomGraph(network)
      rndResult <- callExtremalOptimization(rndNetwork,refine=refine,tau=tau,alpha_max=alpha_max,steps)
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

pcseoss <- function(adjacency,constraints_ml,constraints_cl){
  network <- graph.adjacency(adjacency,mode="undirected",weighted=TRUE)
  res <- callEOPCSEOSS(network,constraints_ml,constraints_cl)
  res[vcount(network)+1] <- round(res[vcount(network)+1],2)
  
  result <- generateOutput(res)

  return(result)
}

###functions for extremal optimization approach
callExtremalOptimization <- function(network,refine,tau,alpha_max,steps){
  
  if(tau==FALSE){
    alpha_max<-1
  }
 
  C_init <- initialWholeGraph(network)
  Q_init <- 0
  
  if(vcount(network)%% 2 == 0){
    half <- vcount(network)/2
  }
  else{
    half <- (vcount(network)+1)/2
  }
  
  sam <-sample(1:vcount(network),half, replace=FALSE)
  G <-NULL
  for(i in 1:vcount(network)){
    if(i %in% sam){
      G[i]=1
    }
    else{
      G[i]=2
    }
  }
  
  C <- getCEO(network,G)
  lambda <- calculateLambda(network,C)
  Q <- calculateQEO(network,lambda)
  
  Q_best <- Q
  C_best <- C
  c_new <- C
  Q_new <- Q
  
  alpha <- 0

  while(alpha<alpha_max){
    #percentage <- alpha/alpha_max*100
    #print(paste(percentage,"%",sep=""))
    
    #generic EO
    vertex <- 0
    if(tau==FALSE){
      vertex <- which.min(lambda)
    }
    #tau-EO
    else if(tau==TRUE){
      lambda_sorted <- sort(lambda, decreasing=FALSE)
      t <- 1+1/log(vcount(network))
      prob <- NULL
      for(i in 1:length(lambda_sorted)){
        prob_i <- i^(-t)
        prob < cbind(prob,prob_i)
      }
      vertex <- sample(1:vcount(network),1,replace=TRUE,prob)
    }

    if(G[vertex]==1){
      G[vertex]=2
    } else if(G[vertex]==2){
      G[vertex]=1
    }
    C_new <- getCEO(network,G)
    #update the fitness function
    lambda <- calculateLambda(network,C_new)
    
    Q_new <- calculateQEO(network,lambda)
    
    if(Q_new>Q_best){
      Q_best <- Q_new
      C_best <- C_new
      alpha <- 0
    } else{
      alpha <- alpha+1
    }
  }

  if(refine=="agent"){
    if(max(C_best)>1){
      tmpRes <- randomLocalSearchAgent(network,C_best,Q_best,steps)
      C_best <- tmpRes[1:vcount(network)]
      Q_best <- tmpRes[vcount(network)+1]
    }
  }
  
  if(Q_best <= 0){
    result <- c(C_init,Q_init)
  }else{
    C <- C_best
    for(c in 1:max(C)){
      vertices <- which(C==c)
      if(length(vertices)>1){
        subgraph <- induced.subgraph(network,vertices)
        alpha_max <- vcount(subgraph)
        tmpResult <- callExtremalOptimization(subgraph,refine=refine,tau=tau,alpha_max=alpha_max,steps=steps)
        C_i <- tmpResult[1:vcount(subgraph)]
        for(i in 1:length(vertices)){
          C_best[vertices[i]]=max(C)+C_i[i]
        }
        C_best <- updateC(C_best)
        Q_best <- calculateQ(network,C_best)
      }
    }
    result <- c(C_best,Q_best)
  }
  return(result)
}

#random local search agents
randomLocalSearchAgent <- function(network,C,Q,maxSteps){
  C_better <- C
  Q_better <- Q
  steps <-0
  while(steps<maxSteps){
    C <- C_better
    c1 <- 0
    c2 <- 0
    v1 <- 0
    v2 <- 0
    while(c1==c2){
      vertices <- sample(1:vcount(network),2,replace=FALSE)
      v1 <- vertices[1]
      v2 <- vertices[2]
      c1 <- C[v1]
      c2 <- C[v2]
    }
    C[v1] <- c2
    C[v2] <- c1
    Q_new <- calculateQ(network,C)
    if(Q_new > Q_better){
      Q_better <- Q_new
      C_better <- C
      steps <- 0
    } else{
      steps <- steps+1
    }
  }
  
  result <- c(C_better,Q_better)
  return(result)
}

#calculate the modularity using the fitness function
calculateQEO <- function(network, lambda){
  Q <- 0
  for(i in 1:vcount(network)){
    Q <- Q+lambda[i]*getDegree(network,i)
  }
  Q <- Q/(sum(get.adjacency(network,attr="weight")))
}

#get community structure for a given bipartition G of the network
getCEO <- function(network,G){
  
  C <- seq(0,0,length.out=vcount(network))
  maxClust <- 0
  
  group_1 <- which(G==1)
  group_2 <- which(G==2)
  
  if(length(group_1>0)){
    sub_1 <- induced.subgraph(network,group_1)
    clusters1 <- clusters(sub_1)$membership
    for(c1 in 1:length(clusters1)){
      vertex <- group_1[c1]
      C[vertex] <- clusters1[c1]
    }
    maxClust <- max(clusters1)
  }
  if(length(group_2>0)){
    sub_2 <- induced.subgraph(network,group_2)
    clusters2 <- clusters(sub_2)$membership
    clusters2 <- sapply(clusters2, function(x) x=x+maxClust)
    for(c2 in 1:length(clusters2)){
      vertex <- group_2[c2]
      C[vertex] <- clusters2[c2]
    }
  }
  return(C)
}

#calculation of the fitness function
calculateLambda <- function(network,C){
  
  lambda <- NULL
  a <- seq(0,0,length.out=max(C))
  
  for(i in 1:vcount(network)){
    comm <- C[i]
    a[comm] <- a[comm]+getDegree(network,i)
  }
  a <- sapply(a, function(x) x=x/(sum(get.adjacency(network,attr="weight"))))
  
  for(i in 1:vcount(network)){
    neighbors <- neighbors(network,i)
    
    numLinks <- 0
    comm <- C[i]
    
    if(length(neighbors)>0){
      for(j in 1:length(neighbors)){
        if(neighbors[j] != i & C[neighbors[j]]==comm){
          numLinks <- numLinks+1
        }
      }
    }

    lambda_v <- numLinks/degree(network,i)-a[comm]
    lambda <- cbind(lambda,lambda_v)
  }
  
  return(lambda)
}

callEOPCSEOSS <- function(network,constraints_ml,constraints_cl){
  #determination of pairwise constraints
  tmp_ml <- rbind(constraints_ml,0)
  tmp_cl <- rbind(constraints_cl,0)
  
  for(m in 1:length(constraints_ml[1,])){
    i <- constraints_ml[1,m]
    j <- constraints_ml[2,m]
    dist <- computeDist(network,i,j)
    tmp_ml[4,m] <- dist
  }
  
  for(c in 1:length(constraints_cl[1,])){
    i <- constraints_cl[1,c]
    j <- constraints_cl[2,c]
    dist <- computeDist(network,i,j)
    tmp_cl[4,c] <- dist
  }
  
  tmp_ml <- tmp_ml[,order(tmp_ml[4,],decreasing=TRUE)]
  
  min_change <- Inf
  threshold <- -1
  
  for(i in 1:length(tmp_ml[4,])){
    thresh <- tmp_ml[4,i]
    ml_change <- length(which(tmp_ml[4,]>thresh))
    cl_change <- length(which(tmp_cl[4,]<=thresh))
    change <- ml_change+cl_change
    if(change<min_change){
      min_change<-change
      threshold <- thresh
    }
  }
  
  constraints_ml <- NULL
  constraints_cl <- NULL

  for(i in 1:length(tmp_ml[1,])){
    if(tmp_ml[4,i]>threshold){
      constraints_cl <- cbind(constraints_cl,tmp_ml[1:3,i])
    }
    else{
      constraints_ml <- cbind(constraints_ml,tmp_ml[1:3,i])
    }
  }
  for(j in 1:length(tmp_cl[1,])){
    if(tmp_cl[4,j]>threshold){
      constraints_cl <- cbind(constraints_cl,tmp_cl[1:3,j])
    }  
    else{
      constraints_ml <- cbind(constraints_ml,tmp_cl[1:3,j])
    }
  }
  
  #extremal optimization using pairwise constraints
  result <- eoPCSEOSS(network,constraints_ml,constraints_cl)
  
  return(result) 
}

eoPCSEOSS <- function(network,constraints_ml,constraints_cl){
  C_init <- initialWholeGraph(network)
  C <- seq(0,0,length.out=vcount(network))
  
  degVertex <- degree(network)
  v_a <- which.max(degVertex)
  v_b <- which.max(degVertex[-v_a])
  if(v_b>=v_a){
    v_b <- v_b+1
  }
  C[v_a] <- 1
  C[v_b] <- 2
  
  for(i in 1:vcount(network)){
    if(i != v_a & i != v_b){
      path_a <- shortest.paths(network,i,v_a,algorithm="dijkstra")
      path_b <-shortest.paths(network,i,v_b,algorithm="dijkstra")
      if(path_a<path_b){
        C[i] <- 1
      }
      else{
        C[i] <- 2
      }
    }
  }

  lambda <- calculateLambdaPCSEOSS(network,C,constraints_ml,constraints_cl)
  Q <- calculateQPCSEOSS(network,C,lambda,constraints_ml,constraints_cl)
  
  Q_best <- Q
  C_best <- C
  alpha <- 1
  
  while(alpha==1){
    alpha <- 0
    
    v <- which.min(lambda)
    if(C[v]==1){
      C[v] <- 2
    }
    else{
      C[v] <- 1
    }
    lambda <- calculateLambdaPCSEOSS(network,C,constraints_ml,constraints_cl)
    Q <- calculateQPCSEOSS(network,C,lambda,constraints_ml,constraints_cl)

    if(Q > Q_best){
      Q_best <- Q
      C_best <- C
      alpha <- 1
    }
  }
  if(Q_best <=0){
    result <- c(C_init,0)
  }
  else{
    C <- C_best
    for(c in 1:max(C)){
      vertices <- which(C==c)
      if(length(vertices)>1){
        subgraph <- induced.subgraph(network,vertices)
        
        #update constraints
        subMl <- constraints_ml
        subCl <- constraints_cl
        
        deleteMl <- NULL
        deleteCl <- NULL
        
        if(length(subMl)>0){
          if(length(subMl)>3){
            for(m in 1:length(subMl[1,])){
              if(!(subMl[1,m]%in%vertices)|!(subMl[2,m]%in%vertices)){
                deleteMl <- c(deleteMl,m)
              }
              else{
                i <- which(vertices==subMl[1,m])
                j <- which(vertices==subMl[2,m])
                subMl[1,m] <- i
                subMl[2,m] <- j
              }
            }
          }
          else{
            if(!(subMl[1]%in%vertices)|!(subMl[2]%in%vertices)){
              subMl <- NULL
            }
            else{
              i <- which(vertices==subMl[1])
              j <- which(vertices==subMl[2])
              subMl[1] <- i
              subMl[2] <- j
            }
          }
          
          if(length(deleteMl)>0){
            subMl <- subMl[,-deleteMl]
          }
        }
        
        if(length(subCl)>0){
          if(length(subCl)>3){
            for(c in 1:length(subCl[1,])){
              if(!(subCl[1,c]%in%vertices)|!(subCl[2,c]%in%vertices)){
                deleteCl <- c(deleteCl,c)
              }
              else{
                i <- which(vertices==subCl[1,c])
                j <- which(vertices==subCl[2,c])
                subCl[1,c] <- i
                subCl[2,c] <- j
              }
            }
          }
          else{
            if(!(subCl[1]%in%vertices)|!(subCl[2]%in%vertices)){
              subCl <- NULL
            }
            else{
              i <- which(vertices==subCl[1])
              j <- which(vertices==subCl[2])
              subCl[1] <- i
              subCl[2] <- j
            }
          }
          
          if(length(deleteCl)>0){
            subCl <- subCl[,-deleteCl]
          }
        }
        
        tmpResult <- eoPCSEOSS(subgraph,subMl,subCl)
        C_i <- tmpResult[1:vcount(subgraph)]
        for(i in 1:length(vertices)){
          C_best[vertices[i]]=max(C)+C_i[i]
        }
        C_best <- updateC(C_best)
        Q_best <- calculateQ(network,C_best)
      }
    }
    result <- c(C_best,Q_best)
  }
  return(result)
}

calculateLambdaPCSEOSS <- function(network,C,constraints_ml,constraints_cl){
  
  lambda <- calculateLambda(network,C)
  
  for(i in 1:length(lambda)){
    if(length(constraints_ml)>3){
      ml <- which(constraints_ml[1,]==i|constraints_ml[2,]==i)
    }
    else{
      ml <- which(constraints_ml[1]==i|constraints_ml[2]==i)
    }
    
    if(length(constraints_cl)>3){
      cl <- which(constraints_cl[1,]==i|constraints_cl[2,]==i)
    }
    else{
      cl <- which(constraints_cl[1]==i|constraints_cl[2]==i)
    }
    
    ml_violation <- 0
    cl_violation <- 0
    
    for(m in ml){
      if(length(constraints_ml)>3){
        i <- constraints_ml[1,m]
        j <- constraints_ml[2,m]
        violation <- constraints_ml[3,m]
      }
      else{
        i <- constraints_ml[1]
        j <- constraints_ml[2]
        violation <- constraints_ml[3]
      }
      
      if(C[i]!=C[j]){
        ml_violation <- ml_violation + violation
      }
    }  
    
    for(c in cl){
      if(length(constraints_cl)>3){
        i <- constraints_cl[1,c]
        j <- constraints_cl[2,c]
        violation <- constraints_cl[3,c]
      }
      else{
        i <- constraints_cl[1]
        j <- constraints_cl[2]
        violation <- constraints_cl[3]
      }
      
      if(C[i]==C[j]){
        cl_violation <- cl_violation + violation
      }
    }
    lambda[i] <- lambda[i]-ml_violation-cl_violation
  }
  
  return(lambda)
}

calculateQPCSEOSS <- function(network,C,lambda,constraints_ml,constraints_cl){
  Q <- calculateQEO(network,lambda)
  
  ml_violation <- 0
  cl_violation <- 0
  
  if(length(constraints_ml)>3){
    for(m in 1:length(constraints_ml[1,])){
      i <- constraints_ml[1,m]
      j <- constraints_ml[2,m]
      if(C[i]!=C[j]){
        ml_violation <- ml_violation + constraints_ml[3,m]
      }
    }
  }
  else{
    if(length(constraints_ml)>0){
      i <- constraints_ml[1]
      j <- constraints_ml[2]
      if(C[i]!=C[j]){
        ml_violation <- constraints_ml[3]
      }
    }
  }
  
  if(length(constraints_cl)>3){
    for(c in 1:length(constraints_cl[1,])){
      i <- constraints_cl[1,c]
      j <- constraints_cl[2,c]
      if(C[i]==C[j]){
        cl_violation <- cl_violation + constraints_cl[3,c]
      }
    }
  }
  else{
    if(length(constraints_cl)>0){
      i <- constraints_cl[1]
      j <- constraints_cl[2]
      if(C[i]==C[j]){
        ml_violation <- constraints_cl[3]
      }
    }
  }
  
  Q <- Q - ml_violation - cl_violation
  return(Q)
}

computeDist<-function(network,i,j){
  distance <- 0
  sum <- 0
  for(v in 1:vcount(network)){
    if(v!=i&v!=j){
      dist_i <- shortest.paths(network,i,v,algorithm="dijkstra")
      dist_j <- shortest.paths(network,j,v,algorithm="dijkstra")
      diff <- (dist_i-dist_j)*(dist_i-dist_j)
      sum <- sum + diff
    }
  }
  sum <- sqrt(sum)
  distance <- sum/(vcount(network)-2)
  
  return(distance)
}
