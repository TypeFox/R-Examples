"diffFree" <-
  function (modellist, diffsfree) 
  {
    ## diffsfree has structure 
    ## list(list(what, ind, dataset, start), ...)
    ## e.g.,  list(list(what="kinpar", ind=c(1,2), 
    ##             dataset=2, start=5), ...)
    ## dataset can be a vector of indices into dataset list
    
    for(diffs in diffsfree){
      for(i in 1:length(diffs$dataset)) {
        if(length(diffs$start) != 0) 
          if(length(diffs$ind) == 2)  
            slot(modellist[[diffs$dataset[i] ]], 
                 diffs$what)[[diffs$ind[1]]][diffs$ind[2]]  <- 
          diffs$start  
        else
          slot(modellist[[diffs$dataset[i] ]], 
               diffs$what)[diffs$ind]  <- diffs$start 
      } 
    }
    modellist
  }

