

mergeMD <- function(list, discard = 0){
  objectlist <- list
  eiform <- objectlist[[1]]$call$formula
  draws.names <- names(objectlist[[1]]$draws)
  for(ii in 2:length(objectlist)){
    if(eiform != objectlist[[ii]]$call$formula | !identical(draws.names, names(objectlist[[ii]]$draws))){stop('eiMD objects not from same model')}
  }


  if(class(objectlist[[1]]$draws[[1]])=="mcmc"){
    samples <- nrow(objectlist[[1]]$draws[[1]])
    for(jj in 1:length(objectlist[[1]]$draws)){
      objectlist[[1]]$draws[[jj]] <-
        objectlist[[1]]$draws[[jj]][(discard + 1):samples,]
      for(ii in 2:length(objectlist)){
        objectlist[[1]]$draws[[jj]] <-
          rbind(objectlist[[1]]$draws[[jj]],
                objectlist[[ii]]$draws[[jj]][(discard + 1):samples,])
      }
      objectlist[[1]]$draws[[jj]] <- mcmc(objectlist[[1]]$draws[[jj]])
    }
  }else{
    samples <- dim(objectlist[[1]]$draws[[1]])[length(dim(objectlist[[1]]$draws[[1]]))]
    for(jj in 1:length(objectlist[[1]]$draws)){
      tmp.names <- dimnames(objectlist[[1]]$draws[[jj]])
      tmp.dims <- dim(objectlist[[1]]$draws[[jj]])
      tmp.dims[length(tmp.dims)] <-
        (tmp.dims[length(tmp.dims)]- discard)*length(objectlist)
      nparam <- prod(tmp.dims[1:(length(tmp.dims)-1)])
      tmp <- as.numeric(objectlist[[1]]$draws[[jj]])[(discard*nparam + 1):(nparam*samples)]
      for(ii in 2:length(objectlist)){
        tmp <- c(tmp, as.numeric(objectlist[[ii]]$draws[[jj]])[(discard*nparam + 1):(nparam*samples)])
      }
      objectlist[[1]]$draws[[jj]] <- array(tmp, c(tmp.dims))
      tmp.names[[length(tmp.names)]] <- 1:tmp.dims[length(tmp.names)]
      dimnames(objectlist[[1]]$draws[[jj]]) <- tmp.names

    }
    
  }


  
  
  objectlist[[1]]$acc.ratios <- NULL
  objectlist[[1]]$sources <- match.call()
  
  return(objectlist[[1]])

  
  
}

