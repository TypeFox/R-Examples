permuteAndScale <- function(Dhat, verbose = FALSE){
  # finds the permutation of Dhat that yields Bhat, having a cycle product smaller than one
  # if such a permutation exsist

  # check whether Dhat already has a CP < 1, if not permute
  if(largestRowElemNotOnDiagonal(Dhat)){
    if(verbose) 
      cat('The largest row element does not lie on diagonal.',
          'Checking whether permutation is stable.\n')
    
    if(possibleToPutLargestRowElemOnDiagonal(Dhat)){
      if(verbose) cat('Permutation is trivial... \n')
      Dhat <- putLargestRowElemOnDiagonal(Dhat)
    }else if(hasCPsmallerOne(Dhat, returnCycleNodes = FALSE)$success){
      if(verbose) 
        cat('Permutation is stable as all cycle products are smaller than one.\n')
    }else{
      if(verbose) cat('Permuting rows to find stable model...\n')
      res.perm <- permute(Dhat, verbose)
      
      if(!res.perm$success){
        warning('Estimate has CP >=1. Model assumptions are not met.', 
                'Returning the empty graph.\n')
      }
      Dhat <- res.perm$Dhat
    }
  }else{
    if(verbose)
      cat('The largest row element lies on diagonal.',
          'So the permutation is stable.\n')
  }
  
  p <- nrow(Dhat)
  
  # scale
  rescaledDhat <- diag(1/diag(Dhat))%*% Dhat
  Ahat <- diag(p) - t(rescaledDhat)
  diag(Ahat) <- 0
  list(Ahat = Ahat, rescaledDhat = rescaledDhat)
}

largestRowElemNotOnDiagonal <- function(Dhat){
  any(diag(abs(Dhat)) - apply(abs(Dhat) - diag(diag(abs(Dhat))), 1, max) <= 0)
}

possibleToPutLargestRowElemOnDiagonal <- function(Dhat){
  max.j.in.rows <- apply(abs(Dhat), 1, which.max)
  !any(duplicated(max.j.in.rows))
}

putLargestRowElemOnDiagonal <- function(Dhat){
  max.j.in.rows <- apply(abs(Dhat), 1, which.max)
  Dhat[order(max.j.in.rows),]
}

hasCPsmallerOne <- function(Dhat, returnCycleNodes, verbose = FALSE){
  # checks whether the cycle product of Bhat is smaller than 1
  # after a few quick checks, we run a variant of Floyd-Warshall 
  # to compute the CP exactly
  
  # with the given Dhat, extract Ahat 
  p <- nrow(Dhat)
  rescaledDhat <- diag(1/diag(Dhat))%*% Dhat
  Ahat <- diag(p) - t(rescaledDhat)
  diag(Ahat) <- 0
  
  # find entries that are >= 1 in absolute value
  entries.larger.than.1 <- as.data.frame(cbind(which(abs(Ahat) >= 1, arr.ind = TRUE), 
                                               coef = abs(Ahat)[abs(Ahat) >= 1]))
  if(verbose) cat('Entries >= 1: ', entries.larger.than.1$coef, '\n')
  
  # ordered list of magnitudes
  ordered.list.of.magnitudes <- sort(abs(Ahat)[abs(Ahat) > 0], decreasing = TRUE)
  
  # if all elements are smaller than 1 or equal to 1, return TRUE
  if(nrow(entries.larger.than.1) == 0){
    if(verbose){
      cat("There are no entries >= 1 in candidate. Thus the model is stable.\n")
    }
    return(list(success = TRUE, cycleNodes = NULL))
  }
  
  # if the elements larger than 1 do not form a dag, 
  # then at least two of them are involved in a cycle
  # with cycle product larger 1 
  # (as we do not consider countering edges with coefficient smaller than 1)
  Ahat.larger.1 <- Ahat
  Ahat.larger.1[which(abs(Ahat.larger.1) < 1)] <- 0
  
  if(!returnCycleNodes){
    G.larger1 <- graph.adjacency(Ahat.larger.1,mode="directed",weighted="a")
    if(!is.dag(G.larger1)){
      if(verbose){
        cat("Edges with weights larger than 1 form a cycle. Model is not stable.\n")
      }
      return(list(success = FALSE, cycleNodes = NULL))
    }
  }
  
  # now, if we add an edge with edge weight smaller than 1
  # it has to be small enough so that a potential cycle
  # does not get a cycle product larger than 1
  
  # the worst case would be that all elements larger than 1 are connected
  # and now adding the largest element with magnitude smaller than 1 closes 
  # the cycle
  # if this worst case does not lead to a cycle product larger than 1,
  # no other case will do so and we have a stable model
  
  product.nodes.larger.one <- abs(prod(Ahat.larger.1[Ahat.larger.1!=0]))
  if(verbose) cat('product', product.nodes.larger.one, '\n')

  if(!returnCycleNodes){
    largest.candidate <- abs(max(ordered.list.of.magnitudes[ordered.list.of.magnitudes < 1]))
    if(largest.candidate*product.nodes.larger.one < 1){
      if(verbose){
        cat("Largest entry smaller than 1", largest.candidate, 
            "will not yield cycle with product >= 1. Model is stable. \n")
      }
      return(list(success = TRUE, cycleNodes = NULL))
    }
  }

  # delete irrelvant nodes:
  small.threshold <- 1/product.nodes.larger.one
  if(verbose) cat('Small threshold:', small.threshold, '\n')
  rows.cols.to.delete <- NULL
  for(i in 1:p){
    max.strength.of.path.though.node.i <- max(abs(Ahat)[i,]) * max(abs(Ahat)[,i])
    if(verbose) cat('Max:', max.strength.of.path.though.node.i, '\n')
    if(max.strength.of.path.though.node.i <= small.threshold){
      rows.cols.to.delete <- c(rows.cols.to.delete, i)
    }
  }

  if(length(rows.cols.to.delete) == 0){
    Ahat.reduced <- Ahat
  }else{
    Ahat.reduced <- Ahat[-rows.cols.to.delete, -rows.cols.to.delete]
  }

  if(verbose){
    cat("Reducing the node set from", p, "to", 
        p-length(rows.cols.to.delete), "nodes and running Floyd-Warshall. \n")
  }
  
  if(p - length(rows.cols.to.delete) <= 1){
    if(verbose){
      cat("All cycle products are smaller than one as only one node', 
          'remained after deleting irrelevant ones.\n")
    }
    return(list(success = TRUE, cycleNodes = NULL))
  }

  # run Floyd-Warshall algorithm and stop as soon as cycle with product >= 1 is found  
  cp.mat <- abs(Ahat.reduced)
  p.relevant <- ncol(Ahat.reduced)
  nodes.in.path.from.i.to.i <- vector("list", length = p.relevant)
  
  for(k in 1:p.relevant){
    for(i in 1:p.relevant){
      for(j in 1:p.relevant){
        path.over.k <- cp.mat[i,k] * cp.mat[k,j]
        path.not.over.k <- cp.mat[i,j]
        
        cp.mat[i,j] <- max(path.not.over.k, path.over.k)
        if(i == j){
          
          if(path.over.k > path.not.over.k){
            nodes.in.path.from.i.to.i[[i]] <- c(nodes.in.path.from.i.to.i[[i]], k)
            if(verbose){
              cat('Going from', i, 'to', j, 'over', k, '.\n')
            }
          }
          
          # if now the cycle product is >= 1, found inadmissible model
          if(cp.mat[i,j] >= 1){
            if(verbose){
              cat("Cycle found with length", cp.mat[i,j], "at node", i, 
                  ". Model is not stable. \n")
              cat('Involved inner nodes are', 
                  nodes.in.path.from.i.to.i[[i]], '.\n')
            }
            nodes.involved.in.cycle <- c(i, nodes.in.path.from.i.to.i[[i]])
            return(list(success = FALSE, cycleNodes = nodes.involved.in.cycle))
          }
        }
      }
    }
  }
  
  if(verbose){
    cat("All cycle products are smaller than one. Largest one is", 
        max(diag(cp.mat)), ".\n", fill = TRUE)
  }
  
  return(list(success = TRUE, cycleNodes = NULL))
}

permute <- function(Dhat, verbose = FALSE){
  # use result from linear assignment problem
  # to find permutation that has a CP as small as possible
  
  p <- nrow(Dhat)
  S <- log(matrix(1,p,p)/abs(Dhat))
  LAP.perm <- as.numeric(solve_LSAP(S))  
  
  Pr <- diag(p)
  Pr <- Pr[,LAP.perm]
  Dhat.tilde <- Pr %*% Dhat
  test <- hasCPsmallerOne(Dhat.tilde, FALSE)
  
  if(test$success){
    if(verbose) cat('The permutation has a cycle product < 1. \n')
    res <- list(Dhat = Dhat.tilde, success = TRUE)
  }else{
    if(verbose) 
      cat('No permutation with cycle product < 1 found.', 
          'Returning the empty graph.\n')
    res <- list(Dhat = diag(p), success = FALSE)
  }
  return(res)
}