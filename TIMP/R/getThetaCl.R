"getThetaCl" <-
  function (th, mod) {
    
    modelspec <- mod@modelspec
    modellist <- mod@modellist
    parorder <- mod@parorder
    datasetind <- mod@datasetind 
    thetaClasslist <- vector("list", length(datasetind))
    for(i in 1:length(modellist))
      thetaClasslist[[i]] <- theta() # consider: thetaClasslist <- lapply( rep("theta", length(modellist)), new )
    for(po in parorder) {
      ds <- po$dataset	
      model <- modellist[[ds[1] ]]
      if(po$name != "prel") {
        parslot <- getPar(model, po, th) 
        p <- po$name
        if(is.list(slot(model, p))) {
          tlist <- list()
          cnt <- 1
          for(i in 1:length(slot(model, p))){ 
            tlist[[i]] <- parslot[cnt:(cnt + 
                                         length(slot(model, p)[[i]]) -1)]
            cnt <- cnt + length(slot(model, p)[[i]])
          }
          parslot <- tlist
        }
        for (d in ds)
          slot(thetaClasslist[[d]], p) <- parslot
      }
    }
    for(i in 1:length(parorder)) {
      if(parorder[[i]]$name == "prel")
        for(j in parorder[[i]]$dataset) {	
          thetaClasslist[[j]] <- addPrelCl(thetaClasslist[[j]], 
                                           modellist[[j]], th, parorder[[i]])
        }
    }
    if(length(mod@data) > 1) 
      thetaClasslist <- getDiffThetaCl(th, thetaClasslist, mod)
    for(i in 1:length(parorder)) {
      if(parorder[[i]]$name == "prel")
        for(j in parorder[[i]]$dataset) {
          thetaClasslist[[j]] <- addPrelCl(thetaClasslist[[j]], 
                                           modellist[[j]], th, parorder[[i]], addM = FALSE)
        }
    }
    thetaClasslist
  }

