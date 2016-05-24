`distance_tm` <-
function(net, projection.method="sum", gconly=TRUE,subsample=1, seed=NULL){
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

  # Project network
  net <- projecting_tm(net, method=projection.method)
  
  # Run the distance calculation, and return output
  return(distance_w(net=net, directed=FALSE,gconly=gconly,subsample=subsample, seed=seed))
}