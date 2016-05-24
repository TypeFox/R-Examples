`rg_reshuffling_tm` <-
function(net,option="links",seed=NULL){
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

  # If seed is set, set it formally
  if(!is.null(seed))
    set.seed(as.integer(seed))

  # Link reshuffling
  if(option=="links") {
    rnet <- cbind(net[,c("i","p")], ok=0)
    E <- nrow(rnet)
    while(sum(rnet[,"ok"])!= E) {
      rE <- which(rnet[,"ok"]==0)
      rnet[rE,2] <- rnet[sample(rE),"p"]
      rnet <- rnet[order(rnet[,"i"], rnet[,"p"]),]
      rnet[,"ok"] <- as.integer(!duplicated(rnet[,c("i","p")]))
      if(sum(rnet[,"ok"])!= E)
        rnet[sample(1:E, size=min(c((E-sum(rnet[,"ok"]))*10, E))),"ok"] <- 0
    }
    rnet <- rnet[,c("i","p")]
    # Add random weights if weighted
    if(ncol(net)==3)
      rnet <- cbind(rnet, w=sample(net[,"w"]))

  # Weight reshuffling
  } else if(option=="weights") {
    if(ncol(net)!=3)
      stop("Weight reshuffling is only possible if the network is weighted")
    rnet <- net
    rnet[,"w"] <- sample(net[,"w"])
  } else {
    stop("Option not recongised, must be either links or weights, see ?rg_reshuffling_tm")
  }
  rownames(rnet)<-NULL;
  attributes(rnet)$tnet <- attributes(net)$tnet
  return(rnet)
}
