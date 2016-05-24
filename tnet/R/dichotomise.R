`dichotomise_w` <-
function(net,GT=0){
  # Ensure that the network conforms to the tnet standard
  if (is.null(attributes(net)$tnet))                      net <- as.tnet(net, type = "weighted one-mode tnet")
  if (attributes(net)$tnet != "weighted one-mode tnet")   stop("Network not loaded properly")

  # Extract ties with a weight greater than GT
  net <- net[net[,"w"]>GT,]
  if(nrow(net)==0) {
    warning("There were no ties with a weight greater than the cutoff")
  } else {
    # Set their weight to 1
    net[,"w"] <- 1
    row.names(net)<-NULL
  }
  return(net)
}

`dichotomise_tm` <-
function(net,GT=0){
  # Ensure that the network conforms to the tnet standard
  if(is.null(attributes(net)$tnet)) {
    if(ncol(net)==3) {
      net <- as.tnet(net, type="weighted two-mode tnet")
    } else {
      warning("This network is a binary two-mode network already")
    }
  }
  if(attributes(net)$tnet!="binary two-mode tnet" & attributes(net)$tnet!="weighted two-mode tnet")
    stop("Network not loaded properly")

  if(attributes(net)$tnet=="weighted two-mode tnet") {
    # Extract ties with a weight greater than GT
    net <- net[net[,"w"]>GT,]
    if(nrow(net)==0) {
      warning("There were no ties with a weight greater than the cutoff")
    } else {
      # Remove weight column
      net <- net[,c("i","p")]
      row.names(net)<-NULL
    }
  }
  attributes(net)$tnet <- "binary two-mode tnet"
  return(net)
}