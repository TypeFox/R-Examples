nodeHarvest <-
function(X,Y, nodesize=10, nodes=1000, maxinter=2, mode="mean", lambda=Inf, addto=NULL,  onlyinter=NULL, silent=FALSE,biascorr=FALSE){
   
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
  
  

  hasnas <- any(is.na(X))
  imputed <- FALSE
  if(hasnas){  
    if(!silent) cat("\n"," imputing missing values for node generation ...")
    tmp <- capture.output( XIMP <- rfImpute(X,if( length(unique(Y))<= 5 ) as.factor(Y) else Y,iter=3)[,-1])
  }
  
  if(!silent) cat("\n ... generating",nodes,"nodes ...")
  Z <- makeRules( if(hasnas) XIMP else X ,Y,nodes=nodes,addZ= addto ,nodesize=nodesize, maxinter=maxinter, onlyinter=onlyinter, silent=silent,levelvec=levelvec)
  tmp <- list(nodes=Z)
  attr(tmp,"levelvec") <- levelvec
  if(hasnas) Z <- adjustmeans(tmp,X,Y)$nodes 
                              
  if(!silent) cat(" ... computing node means ...","\n")
  geti <- getI(Z, if(hasnas) XIMP else X ,Y,mode=mode)
  I <- geti$I
  Z <- geti$Z
  if( any(abs(sapply(Z,function(x) attr(x,"mean")))<10^(-4))){
    geti <- getI(Z, if(hasnas) XIMP else X ,Y=NULL,mode="member")
    Isign <- geti$I
  }else{
    Isign <- abs(sign(I))
  }
  
  wleafs <- rep(0,length(Z))
  indroot <- which(sapply(Z,attr,"depth")==0)[1]
  wleafs[indroot] <- 1
  
  if(!silent) cat(" ... computing node weights ...")
  w <- getw(I,Y,Isign=Isign,wleafs=wleafs, epsilon=lambda-1,silent=silent)
    
  rem <- which(abs(w) < 0.01*max(abs(w)))
  if(length(rem)>0){
    attri <- attributes(Z)
    Z <- Z[-rem]
    attributes(Z) <- attri
    w <- w[-rem]
    I <- I[,-rem,drop=FALSE]
    Isign <- Isign[,-rem,drop=FALSE]
  }
  for (k in 1:length(Z)) attr(Z[[k]],"weight") <- w[k]

  if(ncol(I)>1){
    connection <- t(Isign)%*%Isign
    connection <- diag(1/diag(connection)) %*% connection
    diag(connection) <- 0
    connection[lower.tri(connection)] <- 0
    for (k in 1:nrow(connection)){
      propcontained <- connection[k,]
      maxval <- max(propcontained)
      choose <- which( propcontained>=0.99999 )
      attr(Z[[k]],"ancestors") <- choose
    }
  }else{
    attr(Z[[1]],"ancestors") <- integer(0)
  }
  predicted <- as.numeric(I %*% w)
  
     
  nh <- list()
  nh[["varnames"]] <- colnames(X)
  nh[["nodes"]] <- Z
  nh[["Y"]] <- Y
  class(nh) <- "nodeHarvest"
  if(biascorr){
    corrlin <- lm( Y ~ predicted)
    nh[["predicted"]] <- corrlin$fitted
    nh[["bias"]] <- coef(corrlin)
    if(!silent) cat("  ... applying bias correction (experimental) with coefficients ", coef(corrlin))
  }else{
    nh[["predicted"]] <- predicted
    nh[["bias"]] <- NULL
  }
  return(nh)
   
 }

