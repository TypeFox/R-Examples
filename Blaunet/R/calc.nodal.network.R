calc.nodal.network <-
function(blauObj){

  #initialize
  blauObj$nodalNetwork <- as.data.frame(matrix(0, nrow = nrow(blauObj$memberships), ncol= 2))
  rownames(blauObj$nodalNetwork) <- rownames(blauObj$isInNiche)

  #gets rid of nodes not in the current ecology
  namelist <- network.vertex.names(blauObj$graph)
  diff_names <- setdiff(namelist, rownames(blauObj$dimensions))
  blauObj$graph <- delete.vertices(blauObj$graph, vapply(diff_names, function(x) which(namelist == x), 1))

  edgelist <- as.matrix(blauObj$graph, matrix.type='edgelist')

  #make a named edgelist, makes our computations easier
  charEL <- charEdgelist(edgelist, attr(edgelist, 'vnames'))

  #if we're given an undirected graph (undirected EL/symmetric adjacency matrix)
  #duplicate the EL with the origin nodes reversed
  if (is.directed(blauObj$graph) == FALSE) {
    charEL <- rbind(charEL, cbind(charEL[,2], charEL[,1]))
  }

  #sort edgelist by first element
  if (nrow(charEL) > 1){
    charEL <- charEL[order(charEL[, 1]), ]
  }

  #this is kind of a confusing piece of code at first
  #it sets a 'current' origin node and cycles through all of that node's neighbors
  #when it hits a new 'current' node, it records all of the information for the previous 'current' node
  #then it resets the list of niches spanned to and begins recording information on the new current node
  currentNode <- charEL[1,1]
  spannedTo <- c()

  #cycle through directed edgelist
  #the origin node is element 1, the destination node is element 2
  for (rowCyc in 1:nrow(charEL)){
    edge <- as.vector(charEL[rowCyc,])

    #since EL is sorted, if we see a different origin node,
    #record changes to nodalNetwork
    #update current node
    #reset spannedTo
    if (edge[1] != currentNode){
      blauObj$nodalNetwork[currentNode,1] <- ifelse(length(spannedTo) > 0, 1, 0)
      blauObj$nodalNetwork[currentNode,2] <- length(spannedTo)

      #start new spanner record
      currentNode <- edge[1]
      spannedTo <- c()
      niches1 <- blauObj$isInNiche[edge[1], ]
      niches2 <- blauObj$isInNiche[edge[2], ]
      spannedTo <- union(spannedTo, (which((niches2 - niches1) == 1)))
    }

    else {
      niches1 <- blauObj$isInNiche[edge[1], ]
      niches2 <- blauObj$isInNiche[edge[2], ] 

      #nodal spanners are defined as:
      #node1 is not in nicheA but has a friend in nicheA
      #node1 is then said to 'span' to nicheA

      #niches spanned to are indicated by 1's
      #we get number spanned to
      spannedTo <- union(spannedTo, (which((niches2 - niches1) == 1)))
    }
  }
  
  #save the last elements when loop stops
  blauObj$nodalNetwork[currentNode,1] <- ifelse(length(spannedTo) > 0, 1, 0)
  blauObj$nodalNetwork[currentNode,2] <- length(spannedTo)

  return(blauObj)
}
