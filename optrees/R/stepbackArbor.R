#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Arborescence Problems                                               #
#-----------------------------------------------------------------------------#

# stepbackArbor ---------------------------------------------------------------
#' Go back between two stages of the Edmond's algorithm
#'
#' The \code{stepbackArbor} function rebuilds an arborescence present in 
#' earlier stage of Edmonds's algorithm to find a minimum cost arborescence.
#'
#' @param before list with elements of the previous stage
#' @param after list with elements of the next stage
#'
#' @return A updated list of elements of the earlier stage with a new 
#' arborescence
#' 
#' @seealso This function is an auxiliar function used in 
#' \link{msArborEdmonds} and \link{getMinimumArborescence}.

stepbackArbor <- function(before, after) {
  
  # Previous
  aftStage <- after  # stage with arborescence
  arborArcs <- matrix(aftStage$arborArcs, ncol = 3); arborArcs  # arborescence
  befStage <- before  # previous stage
  superNode <- befStage$superNode; superNode  # supernode that contains cycle
  matchesNodes <- befStage$matches; matchesNodes  # correspondences between nodes
  befArcs <- befStage$arcs; befArcs  # arcs of the previous stage
  befCheapArcs <- befStage$cheapArcs; befCheapArcs
  
  # Four types of arcs
  
  # 1 - Arcs in the new graph already in the previous
  # Arcs not reach or leave the supernode
  i <- which(arborArcs[, 1] != superNode & arborArcs[,2] != superNode)
  arcsType1 <- matrix(arborArcs[i, ], ncol = 3); arcsType1
  
  if (nrow(arcsType1) > 0) {
    # If any arcs found correspondence with arcs of the previous graph
    iniNodes <- matchesNodes[matchesNodes[, 2] %in% arcsType1[, 1], 1]; iniNodes
    endNodes <- matchesNodes[matchesNodes[, 2] %in% arcsType1[, 2], 1]; endNodes
    # Select arcs from previous stage that matches
    newArcs1 <- matrix(ncol = 3)[-1, ]
    for (i in 1:length(iniNodes)) {
      for (j in 1:length(endNodes)) {
        k <- which(befArcs[, 1] == iniNodes[i] & befArcs[,2] == endNodes[j])
        newArcs1 <- rbind(newArcs1, befArcs[k, ])
      }
    }; newArcs1
    
    # If there is more than one arc select the one with minimum cost
    if (any(duplicated(newArcs1[, 2]))) {
      # If there is duplicated nodes
      repNodes <- newArcs1[duplicated(newArcs1[, 2]), 2]; repNodes
      # Remove other rows   
      for (i in 1:length(repNodes)) {
        if (any(newArcs1[, 2] == repNodes[i])) {
          minArc <- min(newArcs1[which(newArcs1[, 2] == repNodes[i]), 3]); minArc
          k <- which(!newArcs1[, 3] == minArc & newArcs1[, 2] == repNodes[i]); k
          if (length(k) == 0) {
            equalArcs <- which(newArcs1[, 2] == repNodes[i]); equalArcs
            newArcs1 <- matrix(newArcs1[-equalArcs[length(equalArcs)], ], ncol = 3); newArcs1
          } else {
            newArcs1 <- matrix(newArcs1[-k, ], ncol = 3); newArcs1
          }
        }
      }
    }
    
    # Add arcs to the arborescence
    befStage$arborArcs <- rbind(befStage$arborArcs, newArcs1)
    
  }
  
  # 2 - Arcs entering the supernode
  i <- which(arborArcs[, 2] == superNode)
  arcsType2 <- matrix(arborArcs[i, ], ncol = 3); arcsType2
  
  if (nrow(arcsType2) > 0) {
    # If there is any arcs we have to find correspondence between nodes
    iniNodes <- matchesNodes[matchesNodes[, 2] %in% arcsType2[, 1], 1]; iniNodes
    endNodes <- matchesNodes[matchesNodes[, 2] %in% arcsType2[, 2], 1]; endNodes
    # Arcs of the previous stages that leaving iniNodes and reaching endNodes
    newArcs2 <- matrix(ncol = 3)[-1, ]
    for (i in 1:length(iniNodes)) {
      for (j in 1:length(endNodes)) {
        k <- which(befArcs[, 1] == iniNodes[i] & befArcs[, 2] == endNodes[j])
        #newArcs2 <- rbind(newArcs2, befArcs[k, ])
        newArcs2 <- rbind(newArcs2, befCheapArcs[k, ])
      }
    }; newArcs2
    
    k <- which(newArcs2[, 3] == min(newArcs2[, 3]))
    mcArc2 <- matrix(newArcs2[k, ], ncol = 3)  # save minimum cost arc
    mcArc2 <- matrix(mcArc2[1,], ncol = 3)  # select first as matrix
    
    # Add arcs to the arborescence
    befStage$arborArcs <- rbind(befStage$arborArcs, mcArc2)
    
  }
  
  # 3 - Arcs leaving supernode
  i <- which(arborArcs[, 1] == superNode)
  arcsType3 <- matrix(arborArcs[i, ], ncol = 3); arcsType3
  
  if (nrow(arcsType3) > 0) {
    # If there is any arcs we have to find correspondence between nodes
    iniNodes <- matchesNodes[matchesNodes[, 2] %in% arcsType3[, 1], 1]
    endNodes <- matchesNodes[matchesNodes[, 2] %in% arcsType3[, 2], 1]
    # Arcs of the previous stages that leaving iniNodes and reaching endNodes
    newArcs3 <- matrix(ncol = 3)[-1, ]
    for (i in 1:length(iniNodes)) {
      for (j in 1:length(endNodes)) {
        k <- which(befArcs[, 1] == iniNodes[i] & befArcs[, 2] == endNodes[j])
        newArcs3 <- rbind(newArcs3, befArcs[k, ])
      }
    }
    
    # If there is more than one arc select the one with minimum cost
    if (any(duplicated(newArcs3[, 2]))) {
      # If there is duplicated nodes
      repNodes <- newArcs3[duplicated(newArcs3[, 2]), 2]  # duplicated nodes
      # Remove other rows
      for (i in 1:length(repNodes)) {
        minArc <- min(newArcs3[which(newArcs3[, 2] == repNodes[i]), 3])
        k <- which(!newArcs3[, 3] == minArc & newArcs3[, 2] == repNodes[i])
        if (length(k) == 0) {
          equalArcs <- which(newArcs3[, 2] == repNodes[i])
          newArcs3 <- matrix(newArcs3[-equalArcs[length(equalArcs)], ], ncol = 3)
        } else {
          newArcs3 <- matrix(newArcs3[-k, ], ncol = 3)
        }
      }
    }
    
    # Add arcs to the arborescence
    befStage$arborArcs <- rbind(befStage$arborArcs, newArcs3)
    
  }
  
  # 4 - Cycle arcs contained within the supernode
  arcsType4 <- befStage$cycleArcs; arcsType4
  
  if (nrow(arcsType4) > 0) {
    # There will be a node already in the arborescence as an arc receptor
    cycleNodes <- matchesNodes[which(matchesNodes[, 2] == superNode), 1]
    receptor <- befStage$arborArcs[befStage$arborArcs[, 2] %in% cycleNodes, 2]
    # Recover the arcs of the cycle
    newArcs4 <- befStage$cycleArcs
    # Remove arc entering receptor node
    newArcs4 <- matrix(newArcs4[-which(newArcs4[, 2] == receptor), ], ncol = 3)
    newArcs4
    # Add arcs to the arborescence
    befStage$arborArcs <- rbind(befStage$arborArcs, newArcs4)
    
  }
  
  # Order the list of arcs of the arborescence
  befStage$arborArcs <- befStage$arborArcs[order(befStage$arborArcs[, 1], 
                                                 befStage$arborArcs[, 2]), ]
  # Recover initial costs of the arcs
  for (i in 1:nrow(befStage$arborArcs)) {
    k <- which(befStage$arborArcs[i, 1] == befArcs[, 1] 
               & befStage$arborArcs[i, 2] == befArcs[, 2])
    befStage$arborArcs[i, 3] <- befArcs[k, 3]
  }; befStage$arborArcs
  
  # Returns list of elements of the stage with updated arborescence
  befStage$arborescence <- TRUE
  befStage$arborArcs
  return(befStage)
  
}
#-----------------------------------------------------------------------------#