weightedID <- function(g, dsc=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dsc))
    dsc <- distSumConnectMatrix(g)

  nn <- numNodes(g)

  cur <- diag(1, nn, nn)
  W_star <- cur
  for (i in 1:(nn - 1)) {
    cur <- cur %*% dsc
    W_star <- W_star + cur
  }

  ID_star <- sum(W_star)
  WID <- nn - 1/nn + ID_star/nn^2

  SID <- sum(diag(W_star))

  list(`WID` = WID, `SID` = SID)
}
