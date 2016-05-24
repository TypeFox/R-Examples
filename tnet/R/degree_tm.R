`degree_tm` <-
function(net, measure = c("degree", "output")){
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

  # Add w=1 column if binary
  if(attributes(net)$tnet=="binary two-mode tnet")
    net <- data.frame(net, w=1)

  # Compute measures
  net <- net[order(net[,"i"], net[,"p"]),]
  out <- data.frame(node=unique(net[,"i"]), degree=NaN, output=NaN)
  index <- cumsum(!duplicated(net[,"i"]))
  if("degree" %in% measure)
    out[, "degree"] <- tapply(net[,"w"], index, length)
  if("output" %in% measure)
    out[, "output"] <- tapply(net[,"w"], index, sum)

  # Add isolates
 if(max(net[,"i"]) != nrow(out)) {
    out <- rbind(out, data.frame(node=1:max(net[,"i"]), degree=0, output=0))
    out <- out[order(out[,"node"]),]
    out <- out[!duplicated(out[, "node"]),]
  }

  # Return table with node ids and chosen measures
  return(out[, c("node", measure)])
}
