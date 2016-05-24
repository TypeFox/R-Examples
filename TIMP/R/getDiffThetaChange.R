"getDiffThetaChange" <-
  function (th, mod) 
  {
    modellist <- mod@modellist 
    parorder <- mod@parorderchange 
    diffchange <- mod@modeldiffs$change   
    for(diffs in diffchange) {
      d <- diffs$dataset[1]
      if(diffs$what %in% slotNames(theta())) {
        fixedpar <- modellist[[d]]@fvecind[[
          which(attributes(modellist[[d]]@fvecind)$names 
                == diffs$what)]]
        ppar <- modellist[[d]]@pvecind[[
          which(attributes(modellist[[d]]@pvecind)$names 
                == diffs$what)]]
        removepar <- sort(unique(append(fixedpar, ppar)))
        parapp <- if(length(removepar) != 0)
          unlist(slot(modellist[[d]], 
                      diffs$what))[-removepar]
        else 
          unlist(slot(modellist[[d]], 
                      diffs$what))
        if(length(parapp) != 0) 
          ind <- (length(th) + 1):(length(th) + length(parapp))
        else ind <- vector()
        parorder[[length(parorder)+1]] <- list(name=diffs$what, 
                                               ind=ind, rm=removepar,  dataset=diffs$dataset)
        
        if(diffs$what %in% modellist[[d]]@positivepar)
          parapp <- log(parapp) 
        else {
          if(length(modellist[[d]]@clinde[[diffs$what]]) > 0)
            for(i in 1:length(modellist[[d]]@clinde[[diffs$what]]))  
              parapp <- log(parapp[i])
          if(length(modellist[[d]]@chinde[[diffs$what]]) > 0) 
            for(i in 1:length(modellist[[d]]@chinde[[diffs$what]]))  
              parapp <- log(parapp[i])
          
        }
        th <- append(th, parapp)
      }
      if(diffs$what == "prelspec") {
        thA <- getPrelTheta(th, modellist, diffs, d, parorder)
        th <- thA$th
        parorder <- thA$parorder
      }
    }
    list(th=th, parorder=parorder)
  }

