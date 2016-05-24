"getXsuper" <-
  function (Xlist, m, t, group)
  {   
    rowLim <- unlist(lapply(Xlist, nrow))
    colLim <- unlist(lapply(Xlist, ncol))
    superX <- matrix(0, sum(rowLim), sum(colLim))
    ## fill in a superX matrix
    for(i in 1:length(group) ){
      rLim <- if(i == 1) 0 else sum(rowLim[1:(i-1)])
      cLim <- if(i == 1) 0 else sum(colLim[1:(i-1)])
      superX[(rLim+1):(rLim+rowLim[i]), 
             (cLim+1):(cLim+colLim[i]) ] <- Xlist[[i]]
      colnames(superX[(rLim+1):(rLim+rowLim[i]), 
                      (cLim+1):(cLim+colLim[i]) ]) <- colnames(Xlist[[i]])
    }
    ## the below adds the constraints to equality based on column labels
    usednames <- vector()
    cnm <- colnames(superX)
    for(i in 1:length(cnm)) {
      nm <- cnm[i]
      if(! nm %in% usednames) {
        usednames <- append(usednames, nm)
        nm_c <- which(colnames(superX) == nm)
        superX[, i] <- apply(superX[, nm_c], 1, sum)
        superX <- superX[, - nm_c[-1] ]
      }
    }
    superX
  }

