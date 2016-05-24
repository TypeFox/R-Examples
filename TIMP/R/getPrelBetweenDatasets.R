"getPrelBetweenDatasets" <- function (modellist, diffsrel) 
{
  ## rel has structure 
  ## list(list(what1, ind1, dataset1,
  ##           what2, ind2, dataset2,
  ##           rel, start, fixed), ...) 
  ## dataset1 can be a vector of indices into dataset list
  ## dataset2 can be a vector of indices into dataset list
  ## length(dataset1) == length(dataset2)
  
  for(diffs in diffsrel){
    for(i in diffs$dataset1) {
      pinde <- slot(modellist[[i]], "mvecind")
      if(length(diffs$ind1) == 1) 
        pinde[[diffs$what1]] <- append(pinde[[diffs$what1]], 
                                       diffs$ind1)
      else 
        pinde[[diffs$what1]] <-
        append(pinde[[diffs$what1]],
               ifelse(diffs$ind1[1] > 1, 
                      length(unlist(slot(modellist[[i]],
                                         diffs$what1)[[1:(diffs$ind1[1] 
                                                          - 1)]]) + diffs$ind1[2]), diffs$ind1[2]))
      slot(modellist[[i]], "mvecind") <- pinde
      if(diffs$what1 == "prel"){
        ind1 <- diffs$ind1
        pspec_ind <- ceiling(ind1/2)
        if(ind1 %% 2 == 0)
          i1 <- ind1 - 1
        else	 i1 <- ind1
        sp <- modellist[[i]]@prelspec[[pspec_ind]]
        modellist[[i]]@nvecind[[sp$what1]] <- append(modellist[[i]]@nvecind[[sp$what1]],  sp$ind1)
        
      }
      
    }
  } 
  modellist
}

