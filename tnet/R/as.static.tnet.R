`as.static.tnet` <-
function(ld){
  # Check that network conforms to tnet standard
  if(is.null(attributes(net)$tnet))
    net <- as.tnet(net, type="longitudinal tnet")
  if(attributes(net)$tnet!="longitudinal tnet")
    stop("Network not loaded properly")

  # Remove time column
  ld <- ld[,c("i","j","w")]

  # Remove joining and leaving of nodes
  ld <- ld[ld[,"i"]!=ld[,"j"],]
  
  # Sort
  ld <- ld[order(ld[,"i"], ld[,"j"]),]

  # Edge index
  index <- !duplicated(ld[,1:2])

  # Create edgelist
  net <- data.frame(ld[index,c("i","j")], w=0)
  dimnames(net)[[2]]<-c("i","j","w")

  # Find weights of ties
  net[,"w"] <- tapply(ld[,"w"], cumsum(index), sum)
  
  # Remove w<=0
  net <- net[net[,"w"]>0,]
  row.names(net)<-NULL
  
  # Check that network conforms to tnet standard
  net <- as.tnet(net,type="weighted one-mode tnet")
  
  # Return object
  return(net)
}