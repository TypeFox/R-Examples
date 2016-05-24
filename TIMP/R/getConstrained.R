"getConstrained" <-
  function (model) 
  {
    # model@constrainedpar is a list of lists 
    # eg
    # contrainedpar = list(list(what = "kinpar", ind = c(1), 
    # low = .3, high = NA))
    
    constrained <- model@constrained 
    ppars <- sort(intersect(slotNames(theta()), slotNames(model)))  
    cinde <- vector("list", length(ppars)) 
    for(f in 1:length(cinde)) cinde[[f]] <- vector()
    names(cinde) <- ppars 
    cindeM <- cinde2M <- cinde2 <- cinde 
    if(length(constrained) > 0) {
      for(r in 1:length(constrained)) {
        if(length(constrained[[r]]$low) > 0) { 
          if(length(constrained[[r]]$ind) == 1)
            c1 <- constrained[[r]]$ind 
          else {
            if(constrained[[r]]$ind[1] > 1) 
              c1 <- (length(unlist(slot(model, 
                                        constrained[[r]]$what)[[1:(
                                          constrained[[r]]$ind[1] - 1)]])) + 
                       constrained[[r]]$ind[2])
            else 
              c1 <- constrained[[r]]$ind[2]
          }
          cinde[[ constrained[[r]]$what ]] <- 
            append( cinde[[ constrained[[r]]$what ]], c1)
          cindeM[[ constrained[[r]]$what ]] <- 
            append( cindeM[[ constrained[[r]]$what ]], 
                    constrained[[r]]$low)
        }
        else ## can't apply constraint from above and below 
          if(length(constrained[[r]]$high) > 0) { 
            if(length(constrained[[r]]$ind) == 1)
              c1 <- constrained[[r]]$ind 
            else {
              if(constrained[[r]]$ind[1] > 1) 
                c1 <- (length(unlist(slot(model, 
                                          constrained[[r]]$what)[[1:(
                                            constrained[[r]]$ind[1] - 1)]])) + 
                         constrained[[r]]$ind[2])
              else 
                c1 <- constrained[[r]]$ind[2]
              
            }
            cinde2[[ constrained[[r]]$what ]] <- 
              append( cinde2[[ constrained[[r]]$what ]], c1)
            cinde2M[[ constrained[[r]]$what ]] <- 
              append( cinde2M[[ constrained[[r]]$what ]], 
                      constrained[[r]]$high)
          }
      }
    }
    model@clinde <- cinde
    model@chinde <- cinde2
    model@lowcon <- cindeM
    model@highcon <- cinde2M
    model
  }

