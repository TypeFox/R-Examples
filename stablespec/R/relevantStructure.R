relevantStructure <- function(listOfFronts, threshold, stableCausal,
                              stableCausal_l1, stableEdge,
                              stringSize) {


  #undirected edges in the inferred causal model
  anyUndirected <- FALSE

  #get at which model complexity the minimum average of BIC is
  minBicAt <- getMinBic(listOfFronts, stringSize)

  #relevant causal paths
  matRelCausal <- stableCausal[[minBicAt]]
  for (i in (minBicAt + 1):length(stableCausal)) {
    matRelCausal <- pmax(matRelCausal, stableCausal[[i]])
  }

  # relevant causal paths length 1
  matRelCausal_l1 <- stableCausal_l1[[minBicAt]]
  for (i in (minBicAt + 1):length(stableCausal_l1)) {
    matRelCausal_l1 <- pmax(matRelCausal_l1, stableCausal_l1[[i]])
  }

  # relevant edges
  matRelEdge <- stableEdge[[minBicAt]]
  for (i in (minBicAt + 1):length(stableEdge)) {
    matRelEdge <- pmax(matRelEdge, stableEdge[[i]])
  }

  #exclude those which lower than threshold
  matRelCausal[which(matRelCausal < threshold)] <-
    matRelCausal_l1[which(matRelCausal_l1 < threshold)] <-
    matRelEdge[which(matRelEdge < threshold)] <- 0

  # for plotting, choose stable edges and direct the edges if appear
  # in causal path length 1
  mat4PlotCausal <- matRelCausal_l1

  mat4PlotCausal[mat4PlotCausal > 0] <- 1

  mat4PlotEdge <- matRelEdge
  mat4PlotEdge[mat4PlotEdge > 0] <- 1

  undirectedEdges <- NULL

  mat4Plot <- mat4PlotEdge + mat4PlotCausal
  for (i in 1:(nrow(mat4Plot) - 1)) {
    #b <- i + 1
    for (j in (i + 1):nrow(mat4Plot)) {

      if (mat4Plot[i, j] != mat4Plot[j, i]) {

        mat4Plot[i, j] <- mat4Plot[i, j] - 1
        mat4Plot[j, i] <- mat4Plot[j, i] - 1

      } else { #insert below line to get the node with undirected edge

          if (mat4Plot[i, j] && mat4Plot[j, i] > 0) {

            undirectedEdges <- rbind(undirectedEdges, c(i, j))
            anyUndirected <- TRUE

          }
      }
    }
  }

  #convert into graph object
  theGraph <- as(mat4Plot, Class="graphNEL")

  #if any undirected edges, then these lines
  #to convert from bi-directed arc to undirected ones
  if (anyUndirected) {
    for (i in 1:nrow(undirectedEdges)) {

      ind <- arrow_h <- NULL
      a <- undirectedEdges[i, 1]
      b <- undirectedEdges[i, 2]
      ind <- paste(a, '~', b, sep="")
      arrow_h[ind] <- "none"
      graph::edgeRenderInfo(theGraph) <-
        list(arrowhead=arrow_h, arrowtail="none")

    }
  }


  # layout the graph
  theGraph <- Rgraphviz::layoutGraph(theGraph)

  return(list(relCausalPath=matRelCausal,
              relCausalPathL1=matRelCausal_l1,
              relEdge=matRelEdge,
              graph=theGraph))
}

getMinBic <- function(listOfFronts, stringSize) {
  # compute the minimum average of BIC
  # To get only a unique model (of a front) in each model complexity
  for (i in 1:length(listOfFronts))
  {
    tempMatA <- listOfFronts[[i]]
    listOfFronts[[i]] <- tempMatA[match(unique(tempMatA[, stringSize + 2]),
                                        tempMatA[, stringSize + 2]), ]
  }

  # Convert list of optimal models into a matrix
  matOfFronts <- do.call(rbind, listOfFronts)

  #order elements in matOfFronts based on their model complexity
  matOfFronts <- matOfFronts[order(matOfFronts[, stringSize + 2]), ]

  #minimum and maximum model complexity
  minCom <- min(matOfFronts[, stringSize + 2])
  maxCom <- max(matOfFronts[, stringSize + 2])
  #theBics <- matOfFronts[, stringSize + 2]
  allBic <- NULL

  for (i in minCom:maxCom) {

    if (is.null(nrow(matOfFronts[which(
      matOfFronts[, stringSize + 2] == i), ]))) {

      #if happens that only a single model complexity exists
      allBic <- c(allBic, matOfFronts[which(
                    matOfFronts[, stringSize + 2] == i), stringSize + 3])
    } else {

      groupedModels <- matOfFronts[which(matOfFronts[, stringSize + 2] == i), ]
      allBic <- c(allBic, mean(groupedModels[, stringSize + 3]))
    }

    # to check whether there is nan values
    # this likely happens if nSubset is extremely small, e.g. 1.
    if (any(is.nan(allBic))) {
      maxBic <- matOfFronts[which(matOfFronts[, stringSize + 2]
                                  == maxCom), stringSize + 3]

      # give nan values the maximum bic of allBic plus 10
      allBic[is.nan(allBic)] <- maxBic + 10
    }
  }
  return(which(allBic == min(allBic)))
}
