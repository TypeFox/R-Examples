`clustering_local_w` <-
function(net, measure="am"){
  # Ensure that the network conforms to the tnet standard
  if (is.null(attributes(net)$tnet))                      net <- as.tnet(net, type = "weighted one-mode tnet")
  if (attributes(net)$tnet != "weighted one-mode tnet")   stop("Network not loaded properly")

  # Find basic parameters
  N <- max(c(net[,"i"],net[,"j"]))
  E <- nrow(net)
  # Ensure network is undirected
  tmp <- symmetrise_w(net, method = "MAX")
  directed <- (nrow(tmp) != nrow(net) | sum(tmp[,"w"]) != sum(net[,"w"]))
  if(directed) 
    stop("Network is not undirected!\nMeasure is not defined from directed networks.\n")
  # Create output object
  net <- net[order(net[,"i"], net[,"j"]),]
  index <- net[,"i"]
  output <- cbind(node=1:N, degree=0, strength=0, am=NaN, gm=NaN, ma=NaN, mi=NaN, bi=NaN)
  output[unique(index), "degree"] <- tapply(net[,"w"], index, length)
  output[unique(index), "strength"] <- tapply(net[,"w"], index, sum)
  # Define numerator-support table
  tri <- cbind(net[,c("i","j")], 1)
  dimnames(tri)[[2]] <- c("j","h","closed")
  # For every node
  for(i in output[output[,"degree"]>=2,"node"]) {
    js <- hs <- net[net[,"i"]==i,c("j","w")]
    dimnames(js)[[2]] <- c("j","wij")
    dimnames(hs)[[2]] <- c("h","wih")
    # All possible ties
    jhs <- merge(js, hs)
    jhs <- jhs[jhs[,"j"]!=jhs[,"h"],]
    jhs <- jhs[,c("j","h","wij","wih")]
    # Find closing ties
    jhs <- merge(jhs, tri, all.x=TRUE)
    jhs[is.na(jhs[,"closed"]),"closed"] <- 0
    jhs <- jhs[,c("wij","wih","closed")]
    jhs <- cbind(jhs, AM=(jhs[,1]+jhs[,2])*0.5,
                      GM=sqrt(jhs[,1]*jhs[,2]),
                      MA=pmax(jhs[,1],jhs[,2]),
                      MI=pmin(jhs[,1],jhs[,2]),
                      BI=1)
    # Calculate ratios
    output[i,"am"] <- sum(jhs[jhs[,"closed"]==1,"AM"])/sum(jhs[,"AM"])
    output[i,"gm"] <- sum(jhs[jhs[,"closed"]==1,"GM"])/sum(jhs[,"GM"])
    output[i,"ma"] <- sum(jhs[jhs[,"closed"]==1,"MA"])/sum(jhs[,"MA"])
    output[i,"mi"] <- sum(jhs[jhs[,"closed"]==1,"MI"])/sum(jhs[,"MI"])
    output[i,"bi"] <- sum(jhs[jhs[,"closed"]==1,"BI"])/sum(jhs[,"BI"])
  }
  # Return output
  return(output[,c("node",measure)])
}



