dissWeights <- function(delta, type = c("unif", "knn", "power", "unifpower"), k = NULL, power = 0){
  # Compute weights as a function of the dissimilarities
  type  <- match.arg(type, c("unif", "knn", "power", "unifpower"), several.ok = FALSE)
  delta <- as.dist(delta)                               # Change delta to a dist class if it is a matrix
  n     <- length(delta)
  n.row <- nrow(as.matrix(delta))
  
  if (type == "unif" | type == "unifpower"){    # such that the weighted empirical distribution (histogram) is uniform
    iord <- order(as.vector(delta),  na.last = TRUE)    # Missing values are ordered last
    y    <- as.vector(delta)[iord]
    y[is.na(y)] <- Inf       # Replace NA by Inf for the computation of tie blocks
    # Find tie blocks
    indTieBlock <- c(1, (2:n)[!y[-n] == y[-1]])         # Index of start of tie block
    ties    <- c(indTieBlock[-1], n + 1) - indTieBlock  # Compute size of tie blocks
    n.ties  <- length(ties)                             # number of different ties
    # Compute midpoints between observed values. 
    mid     <- c(y[indTieBlock[1]],                     # Add smallest y values as first element
                 (y[indTieBlock[-1]] + y[indTieBlock[-length(n.ties)]])/2, 
                 y[indTieBlock[n.ties]])                # Repeat last element 
    width   <- n*(mid[-1] - mid[-n.ties])/(mid[n.ties] - mid[1]) # Compute the width of the interval
    tt      <- rep(width/ties, ties)                    # Down weight width by number of ties and
    # repeat values by number of ties
    w       <- delta                                    # Inherit the class of delta
    w[iord] <- tt                                       # Put values w back in the order of dis
  } 
  
  
  if (type == "knn"){   # Get weights as k nearest neighbors
    w   <- dis.mat <- as.matrix(delta)
    if (is.null(k)) k <- ceil(nrow(w)/4)                # Set the default k to 25% of n
    ind <- apply(dis.mat, 1, FUN = function(delta){sort(delta, index.return = TRUE)$ix[2:(k + 1)]})
    w   <- apply(ind, 2, FUN = function(delta, n.row){w <- rep(0, n.row); w[delta] = 1; return(w)}, n.row)
    w   <- w + t(w)
    w[w > 0] <- 1
    w   <- as.dist(w)
  } 
  
  if (type == "power"){   # Get weights as a power of the dissimilarities
    w   <- delta^power
  }
  
  if (type == "unifpower"){   # First do uniform weighting and then a power of the dissimilarities
    w   <- w * delta^power
  }
  
  return(w)
}
