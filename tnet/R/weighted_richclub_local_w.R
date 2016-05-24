`weighted_richclub_local_w` <-
function(net,prominence){
  # Ensure that the network conforms to the tnet standard
  if (is.null(attributes(net)$tnet))                      net <- as.tnet(net, type = "weighted one-mode tnet")
  if (attributes(net)$tnet != "weighted one-mode tnet")   stop("Network not loaded properly")
  N <- max(c(net[,"i"],net[,"j"]))
  if(length(prominence) != N)
    stop("The length of the prominence parameter is not equal to the number of nodes")
  # Create output table
  output <- cbind(node=1:N, degree=NaN, strength=NaN, degree2p=NaN, num=NaN, den=NaN, ratio=NaN)
  # Go through every node
  for(i in 1:N) {
    # Find the number of ties
    output[i,"degree"] <- sum(net[,1]==i)
    # Find the strength of a node
    output[i,"strength"] <- sum(net[net[,1]==i,3])
    # Find the number of ties to prominent nodes
    output[i,"degree2p"] <- sum(net[,1]==i & net[,2]%in%which(prominence==1))
    # Find the sum of weights to prominent nodes
    output[i,"num"] <- sum(net[net[,1]==i & net[,2]%in%which(prominence==1),3])
    # Close the for-loop for every node
  }

  # Find the randomly expected value
  output[,"den"] <- (output[,"strength"]/output[,"degree"])*output[,"degree2p"]
  # Divide the values on each other to get the ratio
  output[,"ratio"] <- output[,"num"]/output[,"den"]
  # In case no ties to prominent nodes, assign a value of 1
  output[output[,"degree2p"]==0,"ratio"] <- 1
  # Output
  return(output[,c("node","ratio")])
} 
