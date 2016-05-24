mkGLMdf <- function(obj,
                    delta,
                    lwr,
                    upr
                    ) {

  ## check if obj is a list
  if (class(obj)[1] == "list") {
    if (!is.spikeTrain(obj[[1]]) && !is.repeatedTrain(obj[[1]]))
      stop("obj should be a list of spikeTrain or a list of repeatedTrain object.")
    if (is.spikeTrain(obj[[1]])) isST <- TRUE
    else isST <- FALSE
  } else {
    if (!is.spikeTrain(obj) && !is.repeatedTrain(obj))
      stop("obj should be a spikeTrain or a repeatedTrain object.")
    if (is.spikeTrain(obj)) isST <- TRUE
    else isST <- FALSE
    obj <- list(obj)
  } ## End of conditional on class(obj)[1] == "list"

  ## check out lwr
  if (missing(lwr)) {
    if (isST) {
      lwr <- floor(min(sapply(obj,min)))
    } else {
      lwr <- floor(min(sapply(obj,
                              function(l) min(sapply(l,min))
                              )
                       )
                   )
    } ## End of conditional on isST
  } ## End of conditional on missing(lwr)
  ## check out upr
  if (missing(upr)) {
    if (isST) {
      upr <- ceiling(max(sapply(obj,max)))
    } else {
      upr <- ceiling(max(sapply(obj,
                                function(l) max(sapply(l,max))
                                )
                         )
                   )
    } ## End of conditional on isST
  } ## End of conditional on missing(upr)

  ## Find out the number of neurons considered simultaneously
  nbNeurons <- length(obj)
  ## Find out the number of trials
  if (isST) nbTrials <- 1
  else nbTrials <- length(obj[[1]])
  
  ## check out delta
  if (missing(delta)) {
    ## make is smaller than the smallest inter-event interval
    delta <- min(sapply(1:nbNeurons,
                        function(nIdx) {
                          if (isST) {
                            result <- min(diff(obj[[nIdx]]))
                          } else {
                            result <- min(sapply(obj[[nIdx]],
                                                 function(l) min(diff(l))
                                                 )
                                          )
                          }
                          result
                        }
                        )
                 ) - .Machine$double.eps
  }

  ## build the grid
  theGrid <- seq(lwr,upr,delta)
  if (upr > theGrid[length(theGrid)])
    theGrid <- c(theGrid,theGrid[length(theGrid)]+delta)
  theGrid.l <- length(theGrid)
  idxV <- 1:theGrid.l
  
  preNames <- paste("lN.",1:nbNeurons,sep="")
  
  ## define function getTime
  getTime <- function(neuronIdx,trialIdx=NULL) {
    if (isST) as.numeric(obj[[neuronIdx]])
    else as.numeric(obj[[neuronIdx]][[trialIdx]])
  }

  ## define function first list
  firstList <- function(trialIdx=NULL) {

    
    dTM <- matrix(as.integer(0),nrow=theGrid.l,ncol=nbNeurons)
    for (nIdx in 1:nbNeurons) {
      realTime <- getTime(nIdx,trialIdx)
      dTM[findInterval(realTime,theGrid),nIdx] <- 1:length(realTime)
    }
    
    cBRTl <- lapply(1:nbNeurons,
                    function(nIdx) {
                      preMatrix <- matrix(0,
                                          nrow=theGrid.l,
                                          ncol=nbNeurons)

                      for (pIdx in 1:nbNeurons) {
                        evtIdx <- idxV[dTM[,pIdx]>0]
                        if (theGrid.l != evtIdx[length(evtIdx)]){
                          extra <- TRUE
                          evtIdx <- c(evtIdx,theGrid.l)
                        } else {
                          extra <- FALSE
                        }
                        bPts <- diff(evtIdx)
                        ttl <- c(rep(NA,evtIdx[1]),
                                 unlist(lapply(bPts,function(b) 1:b))
                                 )
                        if (!identical(nIdx,pIdx)) {
                          if (extra) ttl[evtIdx[-length(evtIdx)]] <- 0
                          else ttl[evtIdx] <- 0
                        }
                        preMatrix[,pIdx] <- ttl*delta
                      } ## End of for loop on pIdx
                      
                      preMatrix
                      
                    } ## End of function of nIdx
                    )
    dTM[dTM > 0] <- 1
    result <- list(event=as.integer(dTM),
                   time=rep(theGrid,nbNeurons),
                   neuron=factor(rep(1:nbNeurons,each=theGrid.l),
                     levels=1:nbNeurons,
                     labels=paste(1:nbNeurons)
                     )
                   )
    if (!isST) {
      result$trial <- ordered(rep(trialIdx,theGrid.l*nbNeurons),
                              levels=1:nbTrials,
                              labels=paste(1:nbTrials)
                              )
    }

    for (idx in 1:length(preNames)) 
      result[[preNames[idx]]] <- unlist(lapply(cBRTl,function(m) m[,idx]))

    result
  }
  ## End of firstList definition

  if (isST) { 
    result <- as.data.frame(firstList())
  } else {
    result <- lapply(1:nbTrials,firstList)
    varNames <- names(result[[1]])
    result <- lapply(varNames, function(n) unlist(lapply(result, 
            function(l) l[[n]])))
    names(result) <- varNames
    result <- as.data.frame(result)
  } ## End of conditional on isST

  bad <- apply(t(apply(result,1,is.na)),1,any)
  result <- result[!bad,]
  attr(result,"upr") <- upr
  attr(result,"lwr") <- lwr
  attr(result,"delta") <- delta
  attr(result,"call") <- match.call()
  result
  
}
