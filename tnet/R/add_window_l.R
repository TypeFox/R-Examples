`add_window_l` <-
function(net,window=21,remove.nodes=TRUE){
  # Ensure that the network conforms to the tnet standard
  if (is.null(attributes(net)$tnet))                 net <- as.tnet(net, type = "longitudinal tnet")
  if (attributes(net)$tnet != "longitudinal tnet")   stop("Network not loaded properly")
  
  # Set t column as posix
  net[,"t"]       <- as.POSIXct(net[,"t"])
  net <- net[order(net[,"t"]),]
  
  #Check that the fourth column is all 1's
  if(sum(net[,"w"] == rep(1, nrow(net)))!=nrow(net)) 
    stop('The network can only have positive ties (i.e., the fourth column of the data must be 1)');

  #Duplicate all ties, not joining events (i.e., i==j)
  negative       <- net[net[,"i"]!=net[,"j"],];
  
  #Add the window to the timestamps
  negative[,"t"] <- as.POSIXlt(negative[,"t"])+(60*60*24)*window;
  
  #Change the weight column
  negative[,"w"] <- -1
  
  #Merge the two lists
  net <- rbind(net,negative);
  
  #Creating the removal of nodes
  if(remove.nodes) {
    ni <- net[,c("t","i")]
    ni <- ni[nrow(ni):1,]
    ni <- ni[!duplicated(ni[,"i"]),]
    nj <- net[,c("t","j")]
    nj <- nj[nrow(nj):1,]
    nj <- nj[!duplicated(nj[,"j"]),]  
    dimnames(nj)[[2]] <- c("t","i")
    nij <- rbind(ni,nj)
    nij[,"t"] <- as.POSIXlt(nij[,"t"])+(60*60*24)*window;
    nij <- nij[order(nij[,"t"]),]
    nij <- nij[nrow(nij):1,]
    nij <- nij[!duplicated(nij[,"i"]),] 
    nij <- data.frame(nij, j=nij[,"i"], w=-1)
    net  <- rbind(net,nij);
  }
  
  #Order by time
  net <- net[order(net[,"t"], net[,"i"], net[,"j"]),]
  rownames(net) <- NULL;
  return(net)
}