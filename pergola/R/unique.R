#' Creates all trees for distance
#' 
#' Takes order and tree and creates all possible trees
#'  
#' @param hc Hierarchical clustering object.
#' @param dist Distance object.
#' @return Vector with markers with unclear relation (equal distances).
#' 
#' @keywords internal
allTrees <- function(hc, dis, start = 1){
  dis <- as.matrix(dis)
  out <- list(hc$merge)
  for(i in (start + 1):nrow(hc$merge)){
    subl <- getLeaves(hc$merge, hc$merge[i, 1])
    subr <- getLeaves(hc$merge, hc$merge[i, 2])
    pairs <- as.matrix(expand.grid(subl, subr))
    heights <- apply(pairs, 1, function(x) dis[x[1], x[2]])
    sheights <- heights == hc$height[i]
    if(sum(sheights) > 1){
      for(j in which(sheights)){
        hc2 <- switchEntries(hc, e1 = getNode(hc$merge, pairs[j, 1], max = i),
                           e2 = getNode(hc$merge, pairs[j,2], max = i))
        if(nrow(hc$merge) > i){
          out <- c(out, allTrees(hc = hc2, dis = dis, start = i))
        }else{
          out <- c(out, list(hc2$merge))
        }
      }
    }
  }
  return(out)
}

#
#' Get Node for leave
#' 
#' Finds the highest node of a leaf in hclust object below a threshold
#'  
#' @param hc Hierarchical clustering object.
#' @param pair Vector of length 2.
#' @return Updated hclust object.
#' 
#' @keywords internal
getNode <- function(merge, leave, max){
  for(i in max:1){
    leaves <- getLeaves(merge, i)
    if(abs(leave) %in% leaves) return(i)
  }
#   pair <- pair * (-1)
#   hc2$merge[hc$merge == pair[1]] <- pair[2]
#   hc2$merge[hc$merge == pair[2]] <- pair[1]
#   mode(hc2$merge) <- "integer"
#   hc2
  stop("Leave is not in tree.")
}


#
#' Switches entries in hclust object
#' 
#' Switches entries in hclust object
#'  
#' @param hc Hierarchical clustering object.
#' @param pair Vector of length 2.
#' @return Updated hclust object.
#' 
#' @keywords internal
switchEntries <- function(hc, e1, e2){
  if(e1 == e2) return(hc)
  hc2 <- hc
  if(e1 > e2){
    e3 <- e1
    e1 <- e2
    e2 <- e3
  }
  # change smaller entry for merging partner
  arrInd <- which(hc$merge == e1, arr.ind = TRUE)
  e1 <- hc2$merge[arrInd[1], 3-arrInd[1]]
  hc2$merge[hc$merge == e1] <- e2
  hc2$merge[hc$merge == e2] <- e1
  mode(hc2$merge) <- "integer"
  hc2
}



#
#' Subtrees from tree
#' 
#' Provides the subtree from a given tree and node. 
#'  
#' @param merge Tree from hierarchical clustering.
#' @param i Node number.
#' @param leaves Logical. If TRUE the leaves will be provided.
#' @return Vector containing indices.
#' 
#' @keywords internal
getSubtree <- function(merge, i = NULL, leaves = FALSE){
  stack <- merge[i, ]
  if(stack[1] > 0){
    stack <- c(stack, getSubtree(merge, stack[1], leaves))
  }
  if(stack[2] > 0){
    stack <- c(stack,getSubtree(merge, stack[2], leaves))
  }
  if(leaves){
    return(stack)
  }
  stack[stack > 0]
}

#
#' Leaves from subtree
#' 
#' Provides all leaves for a given tree and node.
#'  
#' @param merge Tree from hierarchical clustering.
#' @param i Node number.
#' @return Leave indices.
#' 
#' @keywords internal
getLeaves <- function(merge, i){
  if(i > 0){
    out <- getSubtree(merge, i, TRUE)
    return(abs(out[out < 0]))
  }    
  abs(i)
}

#
#' Calculates SARF 
#' 
#' Internal method to calculate the SARF score for a given order and distance object. 
#'  
#' @param ord Order of markers.
#' @param dis Distance object.
#' @param n Number of neighbors to include.
#' @return SARF score.
#' 
#' @keywords internal
calcSarfDist<-function(ord, dis, n = 1){
  l <- length(ord)
  if(n >= l){
    n <- l - 1
    warning(paste("n was set to ",l - 1))
  }
  m <- as.matrix(dis)[ord, ord]
  sarf <- 0
  for(j in 1:n){
    sarf <- sum(diag(m[, -(1:j)]))
    #for(i in (1 + j):l){    
    #  sarf <- sarf + m[i, i - j]
    #}    
  }
  return(sarf)
}
#
#' Extends SARF criterion to neighborhood
#' 
#' The neighborhood is increased stepwise until only one order remains. 
#'  
#' @param uniOrd Order of markers.
#' @param dis Distance object.
#' @param maxSarf Maximal number of neighbors to include. Should not exceed number of markers -1.
#' @return Order with minimal SARF.
#' 
#' @keywords internal
sarfExt<-function(uniOrd,  dis, maxSarf = length(uniOrd[1, ])-1){
  if(is.vector(uniOrd) || ncol(uniOrd) == 1){
    #warning("Only one order provided.")
    return(uniOrd)
  }
  solved <- FALSE
  for(j in 1:maxSarf){      
    sarfs<-apply(uniOrd, 1, calcSarfDist, dis, j)
    uniOrd<-uniOrd[sarfs == min(sarfs), ]
    if(is.vector(uniOrd)){
      break()
    }
  }
  if(is.vector(uniOrd)){    
    return(uniOrd)
  }else if(nrow(uniOrd) == 2 && all(uniOrd[1, ] == rev(uniOrd[2, ]))){
    if(uniOrd[1, 1] < utils::tail(uniOrd[1, ], n = 1)){
      return(uniOrd[1, ])
    }else{
      return(uniOrd[2, ])
    }
  }  
  else{
    warning(paste("Multiple orders with same SARF score, only top one returned. Increase of maxSarf(", maxSarf, ") might help.", sep = ""))
  }
  return(uniOrd[1, ])
}


