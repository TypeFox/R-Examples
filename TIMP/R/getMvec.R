"getMvec" <-
  function (model) 
  {
    # dependent params. are specified as list or vector
    # in prel; get their indices in the parameter group 
    # in its vectorized form 
    
    ppars <- intersect(slotNames(theta()), slotNames(model))
    pinde <- vector("list", length(ppars))
    for(p in 1:length(pinde)) pinde[[p]] <- vector()
    names(pinde) <- ppars
    
    pinde
  }

