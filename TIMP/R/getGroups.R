"getGroups" <-
  function (m, modeldiffs, labels) 
  {
    thresh <- if (length(modeldiffs$thresh) == 0) 0 else modeldiffs$thresh
    linkclp <- rep(1, length(m))
    sl <- m[[1]]@clpType 
    
    allsl <- datasetind <- clplongind <- vector()
    for (i in 1:length(m)) {
      allsl <- append(allsl, slot(m[[i]], sl))
      datasetind <- append(datasetind, rep(i, length(slot(m[[i]], sl))))
      clplongind <- append(clplongind, 1:length(slot(m[[i]], sl)))
    }
    
    if (length(allsl) < 2) 
      return(list(list(c(1, labels[1]))))
    sort_tmp <- sort(allsl, index.return = TRUE)
    sortclp <- sort_tmp$x
    sortindex <- sort_tmp$ix
    markclp <- sortclp[1]
    groups <- list(list(c(clplongind[sortindex[1]],
                          labels[datasetind[sortindex[1]]])))
    refgroup <- linkclp[datasetind[sortindex[1]]]
    for (i in 2:length(sortclp)) {
      overlimit  <- abs(markclp - sortclp[i]) > thresh 
      if (refgroup != linkclp[datasetind[sortindex[i]]] || overlimit) {
        groups[[length(groups) + 1]] <- list(c(clplongind[sortindex[i]] , 
                                               labels[ datasetind[sortindex[i]]]))
        refgroup <- linkclp[datasetind[sortindex[i]]]
        markclp <- sortclp[i]
      }
      else groups[[length(groups)]][[length(groups[[length(groups)]]) + 
                                       1]] <- c( clplongind[sortindex[i]], labels[datasetind[sortindex[i]] ])
    }
    
    groups
  }
