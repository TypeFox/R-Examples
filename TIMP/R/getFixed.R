"getFixed" <-
  function (model) 
  {
    
    # model@fixed is a list of lists 
    # eg
    # fixed = list(kinpar = list( c(1,2), c(2,3)), irfpar = list(1,2)) 
    # get the fixed parameters indices as vector indices 
    # into the unlisted parameter vector 
    
    fixed <- model@fixed 
    free <- model@free
    ppars <- sort(intersect(slotNames(theta()), slotNames(model)))  
    finde <- vector("list", length(ppars)) 
    for(f in 1:length(finde)) finde[[f]]<-vector()
    names(finde) <- ppars 
    if(length(fixed) > 0) {
      for(r in 1:length(fixed)) {
        if(!is.list(fixed[[r]])){ 
          finde[[ names(fixed)[r] ]] <- as.vector(fixed[[r]]) 
        }
        else {
          for(s in 1:length(fixed[[r]])) {
            indf <- fixed[[r]][[s]][2]
            if(fixed[[r]][[s]][1] > 1){
              slW<-slot(model, names(fixed[r]))
              indf <- indf + sum(unlist(lapply(slW, length))[1:(fixed[[r]][[s]][1]-1)]) 
            }
            finde[[ names(fixed[r]) ]] <-
              append(finde[[ names(fixed[r]) ]], indf)
            
          }
        }
      }
    }
    
    if(length(free) > 0) {
      for(f in 1:length(finde)) finde[[f]]<-vector()
      for(r in 1:length(ppars)) 
        if(length(slot(model, ppars[r])) > 0)
          finde[[ ppars[r] ]] <- 1:length(unlist(slot(model, ppars[r])))
      if(length(free$none) == 0) {
        for(r in 1:length(free)) {
          if(!is.list(free[[r]])){ 
            finde[[ names(free)[r] ]] <- finde[[names(free)[r]]][- 
                                                                   which( finde[[ names(free)[r]]] %in% free[[r]] )]
          }
          else 
            for(s in 1:length(free[[r]])) {
              
              indf <- free[[r]][[s]][2]
              if(free[[r]][[s]][1] > 1){
                slW<-slot(model, names(free[r]))
                indf <- indf + sum(unlist(lapply(slW, length))[1:(free[[r]][[s]][1]-1)]) 
              }
              
              finde[[ names(free)[r] ]] <- finde[[names(free)[r]]][- 
                                                                     which( finde[[ names(free)[r]]] %in% indf)]
            }
        }
      }
    }
    finde
    
  }

