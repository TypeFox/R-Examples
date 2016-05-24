computeBoundary <- function(b1, b0, p, glrTables=NULL, tol=1e-7) {
  p0 <- p[1]
  p1 <- p[2]
  eqSolution <- .computePStar(p, tol)
  N <- ceiling(max(b1, b0) / eqSolution$iStar)
  ##cat("N is ", N, "\n")
  if (is.null(glrTables)) {
    glrTables <- list(g0Table=NULL, g1Table=NULL, tauTable=NULL)
  }
  glrTables <- .updateTables(glrTables, b0=b0, b1=b1, p, N)
  
  g0Table <- glrTables$g0Table
  g1Table <- glrTables$g1Table
  suTable <- rep(NA, N)
  slTable <- rep(NA, N)
  ##cat("B0 is ", b0, "\n")
  for (n in 1:N) {
    ##cat("n is ", n, "\n")
    startS <- ceiling(p0 * n)
    ##cat("start S", startS, "\n")
    searchIndices <- (startS+1):(n+1)
    
    ##cat(g0Table[[n]], "\n")
    k <- which(g0Table[[n]][searchIndices] >= b0)
    ##cat("K is ", k, "\n")
    if ((l <- length(k)) > 0) {
      sSmallest <- searchIndices[k[1]]
      suTable[n] <- sSmallest - 1 ## adjust for 1-based indexing
    }
    endS <- floor(n * p1)
    ##cat("end S", endS, "\n")
    ##cat(g1Table[[n]], "\n")    
    k <- which(g1Table[[n]][1:(endS+1)] >= b1)
    ##cat("K is ", k, "\n")    
    if ((l <- length(k)) > 0) {
      sLargest <- k[l]
      slTable[n] <- sLargest - 1 ## adjust for 1-based indexing
    }
  }

  ##
  ## Trim off the final ends of the table, where suTable == slTable
  ## Looks like the probabilities are zero in any case
  ##
  for (i in ceiling(N/2):N) {
    if (suTable[i] == slTable[i]) break;
  }
  if (i != N) {
    suTable <- suTable[1:i]
    slTable <- slTable[1:i]
    residue <- apply(glrTables$tauTable[(i+1):N, , drop=FALSE], 2, sum)
    if (residue[1] > tol || residue[2] > tol) {
      warning("Boundaries could be inaccurate! Please report to author with complete example. Thanks!")
    }
    N <- i
  }

  ## Re-adjust estimates

  estimate <- apply(glrTables$tauTable[1:N, ], 2, sum)
  names(estimate) <- c("alpha", "beta")

  list(upper=suTable, lower=slTable, estimate=estimate)
}

  
