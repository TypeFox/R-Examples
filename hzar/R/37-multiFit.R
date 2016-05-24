
hzar.multiFitRequest <- function(fitL,each=1,
                                 baseSeed=c(1234,2345,3456,4567,5678,6789,7890,8901,9012,123),
                                 rotateSeed=TRUE,
                                 baseChannel=50,
                                 adjChannel=50,
                                 skip=0){
  
  if(inherits(fitL,"hzar.fitRequest")) {
    fitL <- list(fitL)
  }else if(!is.list(fitL)
           || length(fitL)==0
           || !all(sapply(fitL,inherits,"hzar.fitRequest")))
    stop("fitL is not a list of hzar.fitRequest object.")
  if(!is.numeric(baseChannel))
    baseChannel <- sapply(fitL,function(x) x$mcmcParam$seed[[2]])
  setBaseChannel <- is.numeric(baseChannel)
  setChannelAdj <-  is.numeric(adjChannel) && length(adjChannel)>0 
  
  rotateEverySeed <- (each>1)&&rotateSeed&&!setChannelAdj
  changeEveryChannel <- length(fitL)>1 && (each>1) && (!rotateSeed) && setChannelAdj&& length(baseChannel)==1
  if(changeEveryChannel) setChannelAdj=FALSE
  if(each==1)setChannelAdj=FALSE
  if(!changeEveryChannel&&setBaseChannel)
    baseChannel <- rep(baseChannel,length.out=length(fitL))
  
  if(is.numeric(baseSeed)){
    if(rotateEverySeed){
      baseSeed <- unique(baseSeed)
      L <- length(baseSeed)## baseSeed <- rep(baseSeed,length.out=5+length(fitL)*each)
      seedList <- lapply(skip + 1:(length(fitL)*each),
                         function(x) baseSeed[(cumsum((iter-1) %/% L^(0:5))+0:5) %% L +1])
    }else if(rotateSeed && length(fitL)>1){
      baseSeed <- unique(baseSeed)
      L <- length(baseSeed)## <- rep(baseSeed,length.out=5+length(fitL))
      rotateSeed <- FALSE
      for(iter in 1:length(fitL))
        fitL[[iter]]$mcmcParam$seed[[1]]= baseSeed[(cumsum((skip+iter-1) %/% L^(0:5))+0:5) %% L +1]
    }else{
      baseSeed <- rep(baseSeed,length.out=6)
      rotateSeed <- FALSE
      for(iter in 1:length(fitL))
        fitL[[iter]]$mcmcParam$seed[[1]]= baseSeed[1:6]
          
    }
  }

  if(!changeEveryChannel&&setBaseChannel&&!setChannelAdj){
    for(iter in 1:length(fitL))
      fitL[[iter]]$mcmcParam$seed[[2]]= baseChannel[iter]
  }else if(setChannelAdj&&setBaseChannel){
    channels <- rep(baseChannel[1:length(fitL)],each=each) +
      rep(cumsum(rep(adjChannel,length.out=each)),length(fitL))
    changeEveryChannel=TRUE
  }else if(changeEveryChannel){
    channels <- cumsum(c(baseChannel[1],rep(adjChannel,length.out=each*length(fitL)-1)))
  }
  fitL <- rep(fitL,each=each)
  if(rotateEverySeed)
    for(iter in 1:length(fitL))
      fitL[[iter]]$mcmcParam$seed[[1]]=seedList[[iter]]
  if(changeEveryChannel)
    for(iter in 1:length(fitL))
      fitL[[iter]]$mcmcParam$seed[[2]]=channels[[iter]]
  fitL
  
  
}


  
hzar.doFit.multi <- function(mFitR,doPar=FALSE,inOrder=TRUE){
  fitR=NULL;
  res <- list();
  if(doPar){
    res <- foreach(fitR=mFitR,.inorder=inOrder) %dopar% {
      hzar.doFit(fitR);
    }
  }else {
    res <- foreach(fitR=mFitR,.inorder=inOrder) %do% {
      hzar.doFit(fitR);
    }
  }
  res
}

hzar.doChain.multi <- function(mFitR,doPar=FALSE,inOrder=TRUE,...){
  fitR=NULL;
  res <- list();
  if(doPar){
    res <- foreach(fitR=mFitR,.inorder=inOrder) %dopar% {
      hzar.chain.doSeq(fitR,...);
    }
  }else {
    res <- foreach(fitR=mFitR,.inorder=inOrder) %do% {
      hzar.chain.doSeq(fitR,...);
    }
  }
  res
}
