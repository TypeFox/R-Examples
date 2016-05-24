calc.dyadic <-
function(blauObj, m.dist) { 

  #gets rid of extraneous nodes
  nameList <- network.vertex.names(blauObj$graph)
  diff_names <- setdiff(nameList, rownames(blauObj$dimensions))
  blauObj$graph <- delete.vertices(blauObj$graph, vapply(diff_names, function(x) which(nameList == x), 1))

  edgelist <- as.matrix(blauObj$graph, matrix.type='edgelist')

  charEL <- charEdgelist(edgelist, attr(edgelist, 'vnames'))
  
  #if we're given an undirected graph (undirected EL/symmetric adjacency matrix)
  #duplicate the EL with the origin nodes reversed
  if (is.directed(blauObj$graph) == FALSE) {
    charEL <- unique(rbind(charEL, cbind(charEL[,2], charEL[,1])))
  }

  #sort edgelist by first element
  if (nrow(charEL) > 1){
    charEL <- charEL[order(charEL[, 1]), ]
  }

  if (m.dist == TRUE){
    blauObj$dyadic <- as.data.frame(matrix(0, nrow = nrow(charEL), ncol = 6))
  }
  else{
    blauObj$dyadic <- as.data.frame(matrix(0, nrow = nrow(charEL), ncol = 5))
  }

  edgelistNames <- matrix(0, nrow = 0, ncol = 2)

  #here's where we take advantage of treating the network as directed
  for (rowCyc in 1:nrow(charEL)){
    edge <- as.vector(charEL[rowCyc,])

    edgelistNames <- rbind(edgelistNames, c(edge[1], edge[2]))

	nichea <- blauObj$isInNiche[edge[1],]
	nicheb <- blauObj$isInNiche[edge[2],]
	k <- ncol(blauObj$isInNiche)
	
    if ("ecologyNames" %in% colnames(blauObj$isInNiche)==TRUE) {
	  if (nichea[k]==nicheb[k]) {
	    for (niche in 1:(k-1)) {
		  #CoNicher
		  if (nichea[niche]==nicheb[niche] && nichea[niche]==1) blauObj$dyadic[rowCyc, 1] <- blauObj$dyadic[rowCyc, 1] + 1
		  #spanner
		  if (sum(blauObj$isInNiche[edge[1],(1:(k-1))]) >= 1 && sum(blauObj$isInNiche[edge[2],(1:(k-1))]) >=1 && nichea[niche] + nicheb[niche] == 1) blauObj$dyadic[rowCyc, 4] <- blauObj$dyadic[rowCyc, 4] + 1
		}
		#co-outsider
        if (sum(blauObj$isInNiche[edge[1],(1:(k-1))]) + sum(blauObj$isInNiche[edge[2],(1:(k-1))]) == 0 ){
          blauObj$dyadic[rowCyc, 2] <- 1
        }
		#Straddler
	    if (sum(blauObj$isInNiche[edge[1],(1:(k-1))]) >= 1 && sum(blauObj$isInNiche[edge[2],(1:(k-1))]) == 0){
          blauObj$dyadic[rowCyc, 3] <- sum(blauObj$isInNiche[edge[1],(1:(k-1))])
        }
	    if (sum(blauObj$isInNiche[edge[1],(1:(k-1))]) == 0 && sum(blauObj$isInNiche[edge[2],(1:(k-1))]) >= 1){
          blauObj$dyadic[rowCyc, 3] <- sum(blauObj$isInNiche[edge[2],(1:(k-1))])
        }
	  }
	}  
	else {
	  for (niche in 1:(k)) {
	    #CoNicher
		if (nichea[niche]==nicheb[niche] && nichea[niche]==1) blauObj$dyadic[rowCyc, 1] <- blauObj$dyadic[rowCyc, 1] + 1
		#spanner
		if (sum(blauObj$isInNiche[edge[1],]) >= 1 && sum(blauObj$isInNiche[edge[2],]) >=1 && nichea[niche] + nicheb[niche] ==1) blauObj$dyadic[rowCyc, 4] <- blauObj$dyadic[rowCyc, 4] + 1
	  }	
	  #co-outsider
      if (sum(blauObj$isInNiche[edge[1],]) + sum(blauObj$isInNiche[edge[2],]) == 0 ){
        blauObj$dyadic[rowCyc, 2] <- 1
      }
	  #Straddler
	  if (sum(blauObj$isInNiche[edge[1],]) >= 1 && sum(blauObj$isInNiche[edge[2],]) == 0){
        blauObj$dyadic[rowCyc, 3] <- sum(blauObj$isInNiche[edge[1],])
      }
	  if (sum(blauObj$isInNiche[edge[1],]) == 0 && sum(blauObj$isInNiche[edge[2],]) >= 1){
        blauObj$dyadic[rowCyc, 3] <- sum(blauObj$isInNiche[edge[2],])
      }
	}
	
	#euclidean dist
    blauObj$dyadic[rowCyc,5] <- dist(rbind(blauObj$dimensions[edge[1],], blauObj$dimensions[edge[2],]), method='euclidean')
	
	#mahalanobis dist
    if (m.dist == TRUE){
	  blauObj$dyadic[rowCyc,6] <- sqrt(mahalanobis(blauObj$dimensions[edge[1],], blauObj$dimensions[edge[2],], cov(blauObj$dimensions)))
    }
  }	
  
  blauObj$dyadic <- cbind(edgelistNames, blauObj$dyadic) 
   
  return(blauObj)
}  
	


