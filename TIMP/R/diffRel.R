"diffRel" <-
  function (modellist, diffsrel) 
  {
    ## diffsrel has structure 
    ## list(list(what1, ind1, dataset1,
    ##           what2, ind2, dataset2,
    ##           rel, start, fixed), ...) 
    ## dataset1 can be a vector of indices into dataset list
    ## dataset2 can be a vector of indices into dataset list
    ## length(dataset1) == length(dataset2)
    
    for(diffs in diffsrel){
      for(i in 1:length(diffs$dataset1)) {
        if(length(diffs$rel) == 0 || diffs$rel == "lin"){
          if(length(diffs$ind1)==1 && length(diffs$ind2)==1){
            
            slot(modellist[[diffs$dataset1[i]]], 
                 diffs$what1)[diffs$ind1] <- 
              slot(modellist[[diffs$dataset2[i]]], 
                   diffs$what2)[diffs$ind2] * diffs$start[1] + diffs$start[2]
            
          } 
          if(length(diffs$ind1)==1 && length(diffs$ind2)==2){
            slot(modellist[[diffs$dataset1[i]]], 
                 diffs$what1)[diffs$ind1] <- 
              slot(modellist[[diffs$dataset2[i]]], 
                   diffs$what2)[[diffs$ind2[1]]][diffs$ind2[2]] * 
              diffs$start[1] + diffs$start[2]
            
          }
          if(length(diffs$ind1)==2 && length(diffs$ind2)==1){
            slot(modellist[[diffs$dataset1[i]]], 
                 diffs$what1)[[diffs$ind1[1]]][diffs$ind1[2]] <- 
              slot(modellist[[diffs$dataset2[i]]], 
                   diffs$what2)[diffs$ind2] * diffs$start[1] + diffs$start[2]
            
          }
          if(length(diffs$ind1)==2 && length(diffs$ind2)==2){
            slot(modellist[[diffs$dataset1[i]]], 
                 diffs$what1)[[diffs$ind1[1]]][diffs$ind1[2]] <- 
              slot(modellist[[diffs$dataset2[i]]], 
                   diffs$what2)[[diffs$ind2[1]]][diffs$ind2[2]] * 
              diffs$start[1] + diffs$start[2]
          }
        }
      }       
    }
    modellist
  }

