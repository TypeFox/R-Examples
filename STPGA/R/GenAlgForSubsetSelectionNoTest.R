GenAlgForSubsetSelectionNoTest<-function (P, ntoselect, npop, nelite, mutprob, 
          niterations, lambda, plotiters = TRUE, errorstat="PEVMEAN") 
{
  linenames<-rownames(P)
  InitPop <- lapply(1:npop, function(x) {
    return(sample(linenames, ntoselect))
  })
  if (errorstat=="PEVMEAN"){
   InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
    PEVMEAN(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
  })))
  }
  if (errorstat=="PEVMEAN2"){
    InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
      PEVMEAN2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
    })))
  }
  
  if (errorstat=="PEVMAX"){
    InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
      PEVMAX(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
    })))
  }
  
  if (errorstat=="PEVMAX0"){
    InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
      PEVMAX0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
    })))
  }
  
  if (errorstat=="PEVMAX2"){
    InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
      PEVMAX2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
    })))
  }
  if (errorstat=="CDMEAN"){
    InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
      CDMEAN(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
    })))
  }
  if (errorstat=="CDMEAN2"){
    InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
      CDMEAN2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
    })))
  }
  
  if (errorstat=="CDMAX"){
    InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
      CDMAX(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
    })))
  }
  if (errorstat=="CDMAX0"){
    InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
      CDMAX0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
    })))
  }
  
  if (errorstat=="CDMAX2"){
    InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
      CDMAX2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
    })))
  }
  
  if (errorstat=="DOPT"){
    InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
      DOPT(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
    })))
  }
  
  if (errorstat=="AOPT"){
    InitPopFuncValues <- as.numeric(unlist(lapply(InitPop, FUN = function(x) {
      AOPT(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
    })))
  }
  
  orderofInitPop <- order(InitPopFuncValues, decreasing = FALSE)
  ElitePop <- lapply(orderofInitPop[1:nelite], FUN = function(x) {
    return(InitPop[[x]])
  })
  ElitePopFuncValues <- InitPopFuncValues[orderofInitPop[1:nelite]]
  meanvec <- c()
  for (iters in 1:niterations) {
    CurrentPop <- GenerateCrossesfromElites(Elites = ElitePop, 
                                            Candidates = linenames, npop = npop, mutprob = mutprob)
    if (errorstat=="PEVMEAN"){
    CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                     FUN = function(x) {
                                                       PEVMEAN(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                     })))
    }
    if (errorstat=="PEVMEAN0"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         PEVMEAN0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    if (errorstat=="PEVMEAN2"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         PEVMEAN2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    
    if (errorstat=="PEVMAX"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         PEVMAX(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    if (errorstat=="PEVMAX0"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         PEVMAX0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    
    if (errorstat=="PEVMAX2"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         PEVMAX2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    
    if (errorstat=="CDMEAN0"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMEAN0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    
    if (errorstat=="CDMEAN"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMEAN(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    
    if (errorstat=="CDMEAN2"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMEAN2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    
    if (errorstat=="CDMAX"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMAX(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    
    if (errorstat=="CDMAX0"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMAX0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    
    if (errorstat=="CDMAX2"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMAX2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    
    if (errorstat=="DOPT"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         DOPT(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    
    if (errorstat=="AOPT"){
      CurrentPopFuncValues <- as.numeric(unlist(lapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         AOPT(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda)
                                                       })))
    }
    
    orderofCurrentPop <- order(CurrentPopFuncValues, decreasing = FALSE)
    ElitePop <- lapply(orderofCurrentPop[1:nelite], FUN = function(x) {
      return(CurrentPop[[x]])
    })
    ElitePopFuncValues <- CurrentPopFuncValues[orderofCurrentPop[1:nelite]]
    meanvec <- c(meanvec, min(ElitePopFuncValues))
    if (plotiters) {
      plot(meanvec)
    }
  }
  return(ElitePop)
}


