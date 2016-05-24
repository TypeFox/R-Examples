#the distance matrix needs to have a dim() the same length as the group vector
localoptima <- function(dist, group) {
  c <- length(table(group))
  D <- as.matrix(dist)
  n <- dim(D)[1]
  newgroup <- group
  runs <- 0
  if (any(table(newgroup)==1)) return('There were partitions with a single member. Please check data and try again.')
  while (runs < 1000) {
  #had to insert this for situations with no stable scenarios, ie where the localities all coalesce into 1 group
    if (any(table(newgroup)==1)) return(group)
    #calculate the w/in and w/out group average distances for each sample sample == dim(dist)
    #this needs to be redone every time a locality switches group
    weights<-rclust.weights(newgroup, dist)
    # I now have the weights working out, so now I need to get these things shifting newgroup to find a (local) optimum
    #the following will find all the poorly placed individuals and randomly switch one of them, and repeat the process until all 
    #individuals are optimally placed (find a local optima)
    switchers <- numeric(n)
    for (i in 1:n) {
      if (weights[i,newgroup[i]] != min(weights[i,])) {
        switchers[i] <- 1
      }
    }
    switching <- which(switchers == 1) 
    #if there are any sites that are closer to another group than they are to the one they are currently in
    if (length(switching > 0)) {
      if (length(switching > 1)) {
        switching <- switching[sample(1:length(switching),1)]
        #change these somewhere in here so that it takes only the first if there is more than 1 minimum value
        newgroup[switching] <- which(weights[switching,]==min(weights[switching,]))[1]
      }
      else newgroup[switching] <- which(weights[switching,]==min(weights[switching,]))
      runs <- runs + 1
    }
    else runs <- 1000
  }
  return(newgroup)
}


