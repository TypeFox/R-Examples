####################################################################################
####################################################################################
## AUXILIARY FUNCTIONS TO DEFINE MICROSIMULATION INPUT                            ##
## SZ, November 2013                                                              ##
####################################################################################
####################################################################################

# Set simulation horizon
setSimHorizon <- function(startDate, endDate){
  dts <- c(startDate,endDate)
  simHorizon <- chron(dts,format=c(dates='d/m/Y'),out.format=c(dates='d/m/year'))
  return(simHorizon)
}

# Construct matrix indicating transition pattern and naming the corresponding transition rate functions.
buildTransitionMatrix <- function(allTransitions,absTransitions,stateSpace){  
  if(is.vector(allTransitions))
    allTransitions <- matrix(allTransitions, ncol=2, nrow=1)  
  if(is.vector(absTransitions))
    absTransitions <- matrix(absTransitions, ncol=2, nrow=1)  
  absStates <- absTransitions[,1]  
  if(is.null(dim(stateSpace)))
    stateSpace <- matrix(stateSpace, ncol=1)  
  absStNam <- c('dead')
  if('rest' %in% unlist(strsplit(absStates,'/')))
    absStNam <- c(absStNam, 'rest') 
  transitionMatrix <- matrix(0,nrow=dim(stateSpace)[1], ncol=dim(stateSpace)[1]+length(absStNam))
  colnames(transitionMatrix) <- c(apply(stateSpace,1,paste,collapse='/'),absStNam)
  rownames(transitionMatrix) <- apply(stateSpace,1,paste,collapse='/')  
  # Function to identify whether a set of attributes (`substates') is part of a state space state
  isInThisState <- function(ss,state){
    if(sum(ss %in% as.character(unlist(state)))==length(ss))
      return(TRUE)
    return(FALSE)
  }   
  for(i in 1:length(absStates)){
    strAb <- unlist(strsplit(absStates[i],split='/'))   
    if(length(strAb)==1){
      ia <- which(colnames(transitionMatrix)==absStates[i])
      transitionMatrix[,ia] <- absTransitions[i,2]
    } else {
      iAB <- which(strAb %in% c('dead','rest'))
      aS <- strAb[iAB]
      strAbCov <- strAb[-iAB]
      rA <- which(apply(stateSpace,1,isInThisState, ss=strAbCov)==TRUE)
      ia <- which(colnames(transitionMatrix)==aS) 
      transitionMatrix[rA,ia] <- absTransitions[i,2]
    }
  }        
  if(!is.null(allTransitions)){
    tr <- do.call(rbind,strsplit(allTransitions[,1],'->'))
    for(i in 1: dim(tr)[1]){
      trI <- tr[i,]
      oSPr <- unlist(strsplit(trI[1], split='/'))
      dSPr <- unlist(strsplit(trI[2], split='/'))
      idOS <- apply(stateSpace,1,isInThisState, ss=oSPr)
      idDS <- apply(stateSpace,1,isInThisState, ss=dSPr)
      stateSpaceOS <- stateSpace[idOS,,drop=F]
      stateSpaceDS <- stateSpace[idDS,,drop=F]
      for(j in 1:dim(stateSpaceOS)[1]){
        oS <- as.character(unlist(stateSpaceOS[j,]))
        for(k in 1:dim(stateSpaceDS)[1]){
          dS <- as.character(unlist(stateSpaceDS[k,]))
          c1 <- oS[!oS %in% oSPr] 
          c2 <- dS[!dS %in% dSPr] 
          if(sum(!(c1 %in% c2))==0 & sum(!(c2 %in% c1))==0){
            ir <- which(rownames(transitionMatrix)==paste(oS, collapse='/'))
            ic <- which(colnames(transitionMatrix)==paste(dS, collapse='/'))
            transitionMatrix[ir,ic] <- allTransitions[i,2]
          }          
        }       
      }        
    }   
  }
  return(transitionMatrix)
}