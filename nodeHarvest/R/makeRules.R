makeRules <-
function(X,Y, nodes=1000, addZ=NULL, nodesize=5, maxinter=2, onlyinter=NULL, silent=FALSE, levelvec=NULL){


  if(is.null(levelvec)){

    levelvec <- list()
    for (k in 1:ncol(X)){
      if(!class(X[,k])%in%c("numeric","factor")){
        if(!silent) cat("\n", paste("converting ",class(X[,k]), " variable `",colnames(X)[k],"' to numeric vector in current version ..."),sep="")
        X[,k] <- as.numeric(X[,k])
      }
      if(class(X[,k])=="factor"){
        levelvec[[k]] <- levels(X[,k])
      }else{
        levelvec[[k]] <- numeric(0)
      }
    }
  }
  
  n <- nrow(X)
  ZRULES <- if(!is.null(addZ)) addZ[["nodes"]] else list()
  maxloop <- max(100, round(nodes/10)*2 )
  loopc <- 0
  while(length(ZRULES) < nodes & (loopc <- (loopc + 1))<maxloop){
    rf <- randomForest( X, Y + 0.01*sd(Y)*rnorm(length(Y)), ntree=30, nodesize=nodesize, keep.forest=TRUE, keep.inbag=FALSE, subsample= max(min(round(n/2),10), round(n/5)) , replace=FALSE)
    
    treelist <- list()
    nx <- 0
    
    for (k in 1:rf$ntree){
      tree <- getTree(rf,k=k)
      treelarge <- cbind(tree,rep(0,nrow(tree)))
      colnames(treelarge) <- c(colnames(tree),"level")
      treelarge <- makepostree(treelarge,1,1)
      treelist[[k]] <- treelarge
      nx <- nx + sum(treelarge[,"status"]==-3)
    }
    
    for (treec in 1:length(treelist)){
      tree <- treelist[[ treec ]]
      ZRULES <- c(ZRULES, getRULES(tree,ZRULES,maxinter=maxinter,levelvec=levelvec) )
    }

    
    non <- sapply(lapply(ZRULES,attr,"n"),is.null)
    if(any(non)){
      for (k in which(non)){
        ind <- getsamples(ZRULES[[k]],X,levelvec)
        attr(ZRULES[[k]],"n") <- length(ind)
        attr(ZRULES[[k]],"mean") <- mean(Y[ind]) 
      }
      depth <- sapply(ZRULES,attr,"depth")
      vecn <- sapply(ZRULES,attr,"n")
      
      rem <- numeric(0)
      if( any((vecn < nodesize) & depth>=1) ) rem <- which( (vecn<nodesize) & depth>=1)
      if( any( (vecn > n - nodesize) & depth>=1)) rem <- c(rem, which( (vecn>n-nodesize) & depth>=1))
      rem <- sort(unique(rem))
      if(length(rem)>0){
        ZRULES <- ZRULES[-rem]
      }
    }
    
    if(!is.null(onlyinter)){
      variab <- lapply(ZRULES, function(x) x[,1])
      inter <- sapply(variab,length)
      remove <- !sapply( lapply(variab, function(x) x%in% onlyinter), any) & (inter >1 )
      if(length(remove)>0){
        ZRULES <- ZRULES[-which(remove)]
      }
    }
  }
  attr(ZRULES,"levelvec") <- levelvec

  
  if(length(ZRULES)<0.95* nodes) warning(paste("Could not generate desired ",nodes,"nodes in reasonable time. Working  with",length(ZRULES)," instead"))

  
  if(!any(sapply(ZRULES,attr,"depth")==0)){
    rootnode <- matrix(c(1,-Inf,Inf),nrow=1)
    colnames(rootnode) <- colnames(ZRULES[[1]])
    lZ <- length(ZRULES)
    ZRULES[[lZ+1]] <- rootnode
    attr(ZRULES[[lZ+1]],"n") <- n
    attr(ZRULES[[lZ+1]],"mean") <- mean(Y)
    attr(ZRULES[[lZ+1]],"depth") <- 0
  }
  attri <- attributes(ZRULES)
  ord <- order(sapply(ZRULES,attr,"n"))
  ZRULES <- ZRULES[ord]
  if(!silent) cat("\n", "total number of nodes in initial set                   :", length(ZRULES))
  attributes(ZRULES) <- attri
  
  I <- getI(ZRULES,X,mode="mean")$I
  hash <- round(100000* as.numeric(t(I)%*%rnorm(n)))
  tabb <- table(hash)
  if(any(tabb>1)){
    uniquehash <- as.numeric(names(tabb))
    rem <- numeric(0)
    doub <- which(tabb>1)
    for (khc in 1:length(doub)){
      ind <- which(hash==uniquehash[doub[khc]])
      rem <- c(rem,ind[1:(length(ind)-1)])
    }
    attri <- attributes(ZRULES)
    ZRULES <- ZRULES[-rem]
    attributes(ZRULES) <- attri
  }
  
  if(!silent) cat("\n", "total number of nodes after removal of identical nodes :", length(ZRULES),"\n")
  
  return(ZRULES)
}

