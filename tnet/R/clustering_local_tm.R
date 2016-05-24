`clustering_local_tm` <-
function(net){
  # Ensure that the network conforms to the tnet standard
  if(is.null(attributes(net)$tnet)) {
    if(ncol(net)==3) {
      net <- as.tnet(net, type="weighted two-mode tnet")
    } else {
      net <- as.tnet(net, type="binary two-mode tnet")
    }
  }
  if(attributes(net)$tnet!="binary two-mode tnet" & attributes(net)$tnet!="weighted two-mode tnet")
    stop("Network not loaded properly")

  # Differentiate if binary or weighted two-mode network
  if(attributes(net)$tnet=="binary two-mode tnet")  {
    # Create output table
    output <- data.frame(node=1:max(net[,"i"]), lc=NaN)
    # List of 4-paths
    # 1-paths i1+p1
    paths <- net
    dimnames(paths)[[2]] <- c("i1","p1")
    # 2-paths i1+p1+i2
    dimnames(net)[[2]] <- c("i2","p1")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"i1"] != paths[,"i2"],]
    # 3-paths i1+p1+i2+p2
    dimnames(net)[[2]] <- c("i2","p2")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"p1"] != paths[,"p2"],]
    # 4-paths i1+p1+i2+p2+i3
    dimnames(net)[[2]] <- c("i3","p2")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"i1"] != paths[,"i3"],]
    paths <- paths[paths[,"i2"] != paths[,"i3"],]
    # Find which 4-paths are part of 6-cycles
    dimnames(net)[[2]] <- c("i","p")
    paths <- paths[order(paths[,"i2"], paths[,"i1"], paths[,"i2"], paths[,"p1"], paths[,"p2"]),c("i2","i1","p1","p2","i3")]    
    net.list <- split(net[,"p"], net[,"i"])
    paths <- cbind(paths, closed=apply(paths, 1, function(a) {ct <- c(net.list[[as.character(a["i1"])]],net.list[[as.character(a["i3"])]]); return(sum(duplicated(ct[ct!=as.integer(a["p1"]) & ct!=as.integer(a["p2"])]))>0)}))
    paths <- split(paths[,"closed"], paths[,"i2"])
    for(i in names(paths))
      output[output[,"node"]==i,"lc"] <- sum(paths[[i]])/length(paths[[i]])
  } else {
    # Create output table
    output <- data.frame(node=1:max(net[,"i"]), lc=NaN, lc.am=NaN, lc.gm=NaN, lc.ma=NaN, lc.mi=NaN)
    net <- data.frame(i=as.integer(net[,1]), p=as.integer(net[,2]), w=net[,3])
    # 1-paths i1+p1
    paths <- net
    dimnames(paths)[[2]] <- c("i1","p1","w1")
    # 2-paths i1+p1+i2
    dimnames(net)[[2]] <- c("i2","p1","w2")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"i1"] != paths[,"i2"],]
    # 3-paths i1+p1+i2+p2
    dimnames(net)[[2]] <- c("i2","p2","w3")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"p1"] != paths[,"p2"],]
    # 4-paths i1+p1+i2+p2+i3
    dimnames(net)[[2]] <- c("i3","p2","w4")
    paths <- merge(paths, net, sort=FALSE)
    paths <- paths[paths[,"i1"] != paths[,"i3"],]
    paths <- paths[paths[,"i2"] != paths[,"i3"],]
    pw.am <- as.numeric(rowMeans(paths[,c("w1","w2","w3","w4")]))
    pw.gm <- sqrt(sqrt(paths[,"w1"]*paths[,"w2"]*paths[,"w3"]*paths[,"w4"]))
    pw.ma <- pmax(paths[,"w1"], paths[,"w2"], paths[,"w3"], paths[,"w4"])
    pw.mi <- pmin(paths[,"w1"], paths[,"w2"], paths[,"w3"], paths[,"w4"])
    paths <- paths[,c("i1","p1","i2","p2","i3")]
    # Find which 4-paths are part of 6-cycles
    dimnames(net)[[2]] <- c("i","p","w")
    net <- net[,c("i","p")]
    net.list <- split(net[,"p"], net[,"i"])
    tmp <- cbind(ego=paths[,"i2"], pw.am, pw.gm, pw.ma, pw.mi, closed=apply(paths, 1, function(a) {ct <- c(net.list[[as.character(a["i1"])]],net.list[[as.character(a["i3"])]]); return(sum(duplicated(ct[ct!=as.integer(a["p1"]) & ct!=as.integer(a["p2"])]))>0)}))
    for(i in unique(tmp[,"ego"])) {
      output[output[,"node"]==i,"lc"] <- sum(tmp[,"ego"]==i & tmp[,"closed"])/sum(tmp[,"ego"]==i)
      output[output[,"node"]==i,"lc.am"] <- sum(tmp[tmp[,"ego"]==i & tmp[,"closed"],"pw.am"])/sum(tmp[tmp[,"ego"]==i,"pw.am"])
      output[output[,"node"]==i,"lc.gm"] <- sum(tmp[tmp[,"ego"]==i & tmp[,"closed"],"pw.gm"])/sum(tmp[tmp[,"ego"]==i,"pw.gm"])
      output[output[,"node"]==i,"lc.ma"] <- sum(tmp[tmp[,"ego"]==i & tmp[,"closed"],"pw.ma"])/sum(tmp[tmp[,"ego"]==i,"pw.ma"])
      output[output[,"node"]==i,"lc.mi"] <- sum(tmp[tmp[,"ego"]==i & tmp[,"closed"],"pw.mi"])/sum(tmp[tmp[,"ego"]==i,"pw.mi"])
    }
  }
  # Return output
  return(output)
}



