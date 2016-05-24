labeling <-
function(Y,Q,cd.cluster.object,method = c("2b","2a","1","3"),perm=NULL){

    
    label.method <- match.arg(method)
    #distance <- match.arg(distance)
    
    #------------------------Input check-------------------------#
    input.check(Y,Q,label.method = label.method,perm=perm)
  
    
    Y <- as.matrix(Y)
    Q <- as.matrix(Q)
    
        
    #-----------------------Basic variables----------------------#
    #N:number of examinees
    #J:number of items
    #K:number of attributes
    #M:number of ideal attribute patterns, which is equal to 2^K
    N <- dim(Y)[1]
    J <- dim(Y)[2]
    K <- dim(Q)[2]
    if (label.method=="3"|label.method=="1") {
      M <- 2^K
    }else {
      M <- length(cd.cluster.object$size)
    }
    
    if (M!=length(cd.cluster.object$size)) {
      return(warning("The number of latent clusters cannot be 
                     specified if label method is not '2a' or '2b'."))
    }
    
    #------------------------------------------------------------#
    
  cluster <- as.matrix(cd.cluster.object$class) #cluster memberships
    
    #read alpha
    A <- alpha(K)
    
  Y <- cbind(Y,cluster) #attach cluster membership for each person
  
  matchgroup <- matrix(NA,M,3)
  
    
  #compute the distance matrix between clusters and ideal patterns
  
  #=======================Calculate common weights==============================#
  
  p <- apply(Y[,1:J],2,mean)
  
  if (min(p) == 0 | max(p) == 1)
    {
  
    warning("Cannot compute weights because some weights equal to NA or Inf, unweighted Hamming distance will be used.")
    weight <- c(rep(1,times = J))
  
    }else{
    
      weight <- 1/(p * (1-p))
  
    }
  
  #=======================Calculate distance matrix==============================#
  #dist will be M by M matrix, rows represent clusters and cols represent ideal  #
  #pattern, cell(i,j) represents the distance between ith cluster and jth ideal  #
  #pattern.                                                                      #
  #==============================================================================#
  dist <- NULL
  
  for (i in 1:M)
    {
    
    temp <- DistIdealPatt(Y[which(Y[,(J+1)] == i),1:J],Q,weight)$dist
    dist <- rbind(dist,temp)  #this matrix is ordered from 1 to M
    
    }
  
  #==================================Labeling====================================#
  if (label.method == "2b")
    {
    
    matchgroup[,1] <- c(1:M)    #the number of clusters
    matchgroup[,2] <- as.matrix(apply(dist,1,which.min))#assign labels
  
    } else if (label.method == "1")
    {
    
    # Partially order mean.w based on the order of mean.y
    order <- order(cd.cluster.object$mean.y)
    W <- ((cd.cluster.object$mean.w)[,1:K])[order,]## partial order
    P <- perm
    #change mode
    storage.mode(W) <- "numeric"
    storage.mode(A) <- "integer"
    storage.mode(P) <- "integer"
    K <-as.integer(K)
    length <- as.integer(dim(P)[1])
    #length <- as.integer(2^K)
    res <- matrix(as.integer(0),length,1)
    
    #Compute the min IC
    #dyn.load("ICIndex.dll")
    rownum <- .Fortran("ICIndex",res,K,length,P,A,W)[[1]]
    #Get the min location
    minloc <- which(rownum[,1] == min(rownum),arr.ind = TRUE)
    
    #Get the pattern of min IC
    if (length(minloc) == 1){
      pattern <- P[minloc,]
    }else{
      pattern <- P[minloc[sample(length(minloc),size = 1)],]
    }
    
    matchgroup[,1] <- order
    matchgroup[,2] <- pattern
    
  }else if(label.method == "2a")
  {
    matchgroup[,1] <- c(1:M)
    
    #Compute cluster size
    for (i in 1:M)
    {
      matchgroup[i,3] <-length(Y[which(Y[,(J+1)] == i),(J+1)])
    }
    
    #Order the clusters based on cluster size decreasingly
    #randimize it first in case of ties
    matchgroup <- matchgroup[sample(1:nrow(matchgroup)),]
    matchgroup <- matchgroup[order(matchgroup[,3],decreasing = TRUE,na.last = NA),]

    for (i in matchgroup[,1])
            {
              col <- setdiff(c(1:M),matchgroup[,2])
              if (length(col) != 0)
              {
                if (length(which(dist[i,] == min(dist[i,col]),arr.ind = TRUE)) == 1) #no ties
                {
                  matchgroup[which(matchgroup[,1] == i,arr.ind = TRUE),2] <- which(dist[i,] == min(dist[i,col]),arr.ind = TRUE)
                }else #there exist ties
                {
                  x <- setdiff(which(dist[i,] == min(dist[i,col]),arr.ind = TRUE),matchgroup[,2])
                  if (length(x) == 1){
                    matchgroup[which(matchgroup[,1] == i,arr.ind = TRUE),2] <- x
                  } else {
                    matchgroup[which(matchgroup[,1] == i,arr.ind = TRUE),2] <- sample(x,size = 1)
                  }
                  
                }
              }
            }
    
  } else if (label.method == "3"){
    #Compute Mean for each cluster
    meanScore <- aggregate(apply(Y[,1:J],1,sum),by = list(Y[,(J+1)]),FUN = "mean")
    
    matchgroup[,1] <- meanScore[,1]#number of cluster
    matchgroup[,3] <- meanScore[,2]#size of cluster
    
    #randimize it in case of ties
    matchgroup <- matchgroup[sample(1:nrow(matchgroup)),]
    matchgroup <- matchgroup[order(matchgroup[,3]),]
    
    matchgroup[1,2] <- 1 #min value
    matchgroup[M,2] <- M #max value
    
    #========================================IF nested================================================#
    # test if one cluster is nested within another one cell(i,j)=0 represents ith cluster is NOT nested      
    # within jth cluster. cell(i,j)!=0 represents ith cluster is nested with jth cluster.      
    #=================================================================================================#
    if.nested <- alpha(K)%*%t(alpha(K))   
    
    
    for (i in 2:M){
      for (j in 2:M){
        if (if.nested[i,j] == max(if.nested[i,])){
          if.nested[i,j] <- 1
        }else{
          if.nested[i,j] <- 0
        }
      }
    }
    
    #========================initialize gi========================#
    #       gi is the clusters of which rowsum of alpha is 1      #
    #       gi is the group in which person only master one skill #
    #=============================================================#
    index <- matrix(NA,M,2)
    index[,1] <- as.matrix(seq(1:M))
    index[,2] <- as.matrix(apply(alpha(K),1,sum))
    gi <- index[which(index[,2] == 1),1]
    
    
    for (i in 2:(M-1))
    { 
      #================finalize gi========================#
      gi <- setdiff(gi,matchgroup[,2])
      #look for gplus
      
      for (gplus in 2:M){
        if (apply(as.matrix(if.nested[gi,gplus]),2,sum) == 0){ #gi is not nested within gplus
          gi <- setdiff(union(gi,gplus),matchgroup[,2]) #update gi
        }
      }
      
      #================labelling========================#
      g <- gi #all possible clusters for labelling at current step
      
      if (length(g) != 0)
      {
        if (length(which(dist[matchgroup[i,1],g] == min(dist[matchgroup[i,1],g]),arr.ind = TRUE)) == 1) # no ties
        {
          matchgroup[i,2] <- g[which(dist[matchgroup[i,1],g] == min(dist[matchgroup[i,1],g]),arr.ind = TRUE)]
        }else #there are ties
        {
          x <- which(dist[matchgroup[i,1],g] == min(dist[matchgroup[i,1],g]),arr.ind = TRUE)
          if (length(x) == 1){
            matchgroup[i,2] <- g[x]
          }else{
            matchgroup[i,2] <- g[sample(x,size = 1)]
          }
          
        }
        
      }
    }
    if (length(gi) == 0){ #no group need to be labelled
      break
    }
  }
    #================Calculate output================#
    label.cluster <- array(NA,dim = length(cluster))
    for (i in 1:length(cluster)){
      if (cluster[i]%in%matchgroup[,1]){
        label.cluster[i] <- matchgroup[which(matchgroup[,1] == cluster[i]),2]
      }
    }
    
    label <- matrix(NA,N,K)
    for (i in 1:N)
    {
      label[i,] <- alpha(K)[label.cluster[i],]
    }
    
    att.patt <- rep(NA,times=2^K)
    if (M==2^K){
      freq <- rep(0,times=2^K)
    } else {
      freq <- rep(NA,times=2^K)
    }
    
    for (i in 1:2^K){
      att.patt[i] <- paste(A[i,],collapse=" ")
      freq[i] <- sum(label.cluster==i)
    }
    att.dist <- data.frame(att.patt,freq)
    output <- list(att.pattern = label, att.class = label.cluster, att.dist = att.dist, label.method=label.method)
    class(output) <- "labeling"
  return(output)
}
