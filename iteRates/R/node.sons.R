node.sons <-
function(phy, node)#Thanks to Julien Dutheil - code from earlier version of APE
{
  if (!("phylo" %in% class(phy))) stop("Object \"phy\" is not of class \"phylo\"") 
  E <- phy$edge
  n <- dim(E)[1]
  sons <- numeric(0)
  count <- 1
  for(i in 1:n) {
    if(E[i,1] == node) {
      sons[count] <- E[i,2];
      count <- count + 1
    }
  } 
  return(sons)
}

