netProp = function(network, from = "FROMJCT", to = "TOJCT", pred = "pred", iprint = 1) {
  warning("This function is deprecated and will be removed in a future version of rtop")
  if (requireNamespace("igraph")) {

    rndf = data.frame(FROMJCT = network$FROMJCT, TOJCT = network$TOJCT, 
               OBJECTIT = network$OBJECTID, pred = network@data[,pred])
    igr = igraph::graph.data.frame(rndf)
    igrs = igraph::topological.sort(igr, mode = "out")
    rndf$to = match(as.character(rndf$TOJCT), igraph::V(igr)$name[igrs+1])
    while (TRUE) {
      lcon = which(rndf$to == min(rndf$to[is.na(rndf$pred)]))
      if (length(lcon) == 0) 
        break()
      if (iprint > 0) print(lcon)
      while (is.na(rndf$pred[lcon[1]])) {
        ncon = igraph::neighbors(igr, lcon[1] - 1) + 1
        if (length(ncon) == 0 || ncon > dim(rndf)[1]) {rndf$pred[lcon[1]] = -9999; break}
        lcon = c(ncon, lcon)
        print(lcon)
      }
      rndf$pred[lcon] = rndf$pred[lcon[1]]
    }
    rndf$pred[rndf$pred < -9998] = NA
    network@data[, pred] = rndf$pred
  } else {
    warning("This function will perform faster with igraph installed")
    ichange = 1
    while (ichange > 0) {
      ichange = 0
      for (i in 1:dim(network)[1]) {
        if (!is.na(network@data[i,pred])) {
          tt = which(network@data[,to] == network@data[i,from])
          if (length(tt) > 0 && is.na(network@data[tt,pred])) {
            network@data[tt, pred] = network@data[i, pred]
            ichange = ichange + 1
          }
        }
      }
      if (iprint > 0) print(paste("ichange",ichange))
    }
  }
  network
}
