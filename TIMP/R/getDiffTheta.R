"getDiffTheta" <- function(th, mod) {
  modellist <- mod@modellist 
  modeldiff <- mod@modeldiffs
  parorder <- list()
  if(length(modeldiff$change) != 0) {
    thA <- getDiffThetaChange(th, mod)
    th <- thA$th
    mod@parorderchange <- thA$parorder
  }
  if(length(modeldiff$free)!=0 || length(modeldiff$add) !=0)
    for(diff in append(mod@modeldiffs$free, mod@modeldiffs$add)){ 
      j <- diff$dataset[1] 
      fixed <- modellist[[j]]@fvecind
      prel <- modellist[[j]]@pvecind 
      model <- modellist[[diff$dataset[1]]]
      removepar <- sort(append(fixed[[diff$what]],prel[[diff$what]]))
      if(length(diff$ind) == 2) {
        whichfree <- diff$ind[2]
        if(diff$ind[1] > 1) {
          slW<-slot(model, diff$what)
          whichfree <- whichfree + sum(unlist(lapply(slW, length))[1:(diff$ind[1]-1)]) 
        }
        partmp <- slot(model, diff$what)[[diff$ind[1]]][diff$ind[2]]
        
      }
      else { 
        partmp <- slot(model, 
                       diff$what)[diff$ind]
        whichfree <- diff$ind
      }
      if(diff$what %in% model@positivepar)
        partmp <- log(partmp)
      else {
        if(length(model@clinde$po$name) > 0) 
          for(i in 1:length(model@clinde[[diff$what]]))  
            parapp <- log(parapp[i])
        if(length(model@chinde[[diff$what]]) > 0) 
          for(i in 1:length(model@chinde[[diff$what]]))  
            parapp <- log(parapp[i])
      }
      removepar <- which(whichfree %in% removepar)
      if(length(removepar)>0)
        partmp <- partmp[-removepar]
      if(length(partmp) > 0)
        ind <- (length(th) + 1):(length(th) + length(partmp)) 
      else ind <- vector()  
      parorder[[length(parorder)+1]] <- list(name=diff$what,
                                             ind=ind, rm=removepar, dataset=diff$dataset, indm=diff$ind)
      th <- append(th, partmp)
      
    }          
  if(length(modeldiff$rel)!=0)
    for(diff in modeldiff$rel) {
      removepar <- if(diff$fix) 1:length(diff$start) else vector()
      if(length(removepar)>0)
        partmp <- diff$start[-removepar]
      ind <- if(length(partmp) > 0) (length(th) +
                                       1):(length(th) + length(partmp)) else vector()
      parorder[[length(parorder)+1]] <- list(name=diff$what1,
                                             ind=ind, rm=removepar, dataset=diff$dataset1, indm=diff$ind1)
      th <- append(th, partmp)
    }
  if(length(modeldiff$dscal) != 0) {
    for(i in 1:length(modeldiff$dscal)) {
      j <- modeldiff$dscal[[i]]$to 
      if(length(modellist[[j]]@drel) != 0){ 
        fixed <- modellist[[j]]@fvecind
        prel <- modellist[[j]]@pvecind 
        removepar <- sort(append(fixed[["drel"]], prel[["drel"]]))
        parapp <- if(length(removepar)!=0)  
          unlist(slot(modellist[[j]], 
                      "drel"))[-removepar] 
        else unlist(slot(modellist[[j]], 
                         "drel"))
        if(length(parapp) != 0)
          ind <- (length(th) + 1):(length(th) + length(parapp))
        else ind <- vector()
        parorder[[length(parorder)+1]] <- list(name="drel",
                                               ind=ind, rm=removepar, dataset=j, 
                                               indm=(1:length(modellist[[j]]@drel)))
        th <- append(th, parapp)
      }
    }
  }
  mod@parorderdiff <- parorder
  list(theta = th, mod = mod)
}

