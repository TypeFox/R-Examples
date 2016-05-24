perc.cutoffs <-
function(CnsIndx, OrdIndx, Y, N){
    perc.cut <- list()
    for(j in (CnsIndx+1):OrdIndx){
      perc.cut[[j]] <- qnorm(c(0, cumsum(table(Y[, j])/N)))
    }
    perc.cut
  }
