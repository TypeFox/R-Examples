glrSearch <-
function(p, alpha, beta, stepSize=0.5, tol=1e-7,
                      startB1 = log(1/beta), startB0 = log(1/alpha),
                      maxIter=25, gridIt=FALSE, nGrid=5,
                      verbose=FALSE) {
  if ((length(p) != 2) || (p[1] >= p[2])) {
    stop("glrSearch: The p vector should be c(p0, p1) with p0 < p1!")
  }

  count <- 1

  ## Compute tables and \hat{alpha} and \hat{beta}, estimates of alpha and beta
  glrTables <- list(g0Table=NULL, g1Table=NULL, tauTable=NULL)

  ##
  ## Start searching for b and a
  ##
  ## Solve for pStart
  eqSolution <- .computePStar(p, tol=tol)

  ## Our starting a and values
  b1 <- startB1; b0 <- startB0

  N <- ceiling(max(b1, b0) / eqSolution$iStar)


  glrTables <- .updateTables(glrTables, b0=b0, b1=b1, p, N)
  estimate <- apply(glrTables$tauTable, 2, sum)
  names(estimate) <- c("alpha", "beta")
  resultList <- list(b1=log(1/alpha), b0=log(1/beta), estimate=estimate)

  if (verbose) {
    print(paste("New N", N))
    print(paste("Trying b0:", b0, "and b1:", b1))
  }

  banging <- function(resultList) {
    n <- length(resultList$b1)
    if (n <= 2) {
      FALSE
    } else {
      b1.2 <- resultList$b1[n]; b0.2 <- resultList$b0[n];
      b1.0 <- resultList$b1[n-2]; b0.0 <- resultList$b0[n-2];
      (sqrt((b1.0 - b1.2)^2 + (b0.0 - b0.2)^2) < tol)
    }
  }

  ##
  ## Keep last two a's and b's to prevent oscillation
  ##
  repeat {
    count <- count + 1
    if (count > maxIter) {
      ## We have exhausted our tries
      break;
    }

    absBias  <- abs(estimate - c(alpha, beta))

    if (sqrt(sum(absBias^2)) < tol) {
      ## we are done
      break;
    }

    ## Stop if you are banging against two limits
    if (banging(resultList)) {
        if (verbose) print("Yeah, Banging between two limits")
        break;
    }

    ## Otherwise Do a greedy reduction
    if (absBias[1] > absBias[2]) {
      ## try a new B0
      if (verbose) print("Trying a new B0 to refine alpha")
      b0 <- b0 + ifelse(estimate[1] > alpha, stepSize, -stepSize)
    } else {
      ## try a new B1
      if (verbose) print("Trying a new B1 to refine beta")
      b1 <- b1 + ifelse(estimate[2] > beta, stepSize, -stepSize)
    }
    N <- ceiling(max(b0, b1) / eqSolution$iStar)
    glrTables <- .updateTables(glrTables, b0=b0, b1=b1, p, N)
    estimate <- apply(glrTables$tauTable, 2, sum)
    names(estimate) <- c("alpha", "beta")

    if (verbose) {
      print(paste("Trying b0:", b0, "and b1:", b1))
      print(paste("New N", N))
      print("Estimates:")
      print(estimate)
    }
    resultList <- list(b1=c(resultList$b1, b1), b0=c(resultList$b0, b0),
                       estimate=rbind(resultList$estimate, estimate))
  }

  ##browser()
  ## Note down the count
  n <- length(resultList$b1)

  ##
  ## If banging, we need bracket the other variable
  ##
  ##stop("Testing")

  if (banging(resultList)) {
      b1.1 <- resultList$b1[n-1]; b0.1 <- resultList$b0[n-1];
    estimate1 <- resultList$estimate[n-1, ]

    n2 <- n
    if (b1.1 != b1) {
      b1Limits <- c(b1.1, b1)
      if (verbose) print("b1 values were banging...")
      ## change b0 after picking the b1 that goes below beta
      if (estimate1["beta"] <= beta) {
        b1 <- b1.1
        b0 <- b0.1
        n2 <- n-1
      } # ow b1 and b0 remain the same..
      b0Limits <- c(b0, 0)
      if (estimate1["alpha"] <= alpha) {
        repeat {
          count <- count + 1
          if (count > maxIter) {
            ## We have exhausted our tries
            break;
          }
          b0 <-  b0 - stepSize
          N <- ceiling(max(b0, b1) / eqSolution$iStar)
          glrTables <- .updateTables(glrTables, b0=b0, b1=b1, p, N)
          estimate <- apply(glrTables$tauTable, 2, sum)
          names(estimate) <- c("alpha", "beta")
          if (verbose) {
            print("B0 going down")
            print(paste("Trying b0:", b0, "and b1:", b1))
            print(paste("New N", N))
            print("Estimates:")
            print(estimate)
          }
          resultList <- list(b1=c(resultList$b1, b1), b0=c(resultList$b0, b0),
                             estimate=rbind(resultList$estimate, estimate))
          if (estimate["alpha"] >= alpha) {
            break;
          }
        }
      } else {
        repeat {
          count <- count + 1
          if (count > maxIter) {
            ## We have exhausted our tries
            break;
          }
          b0 <-  b0 + stepSize
          N <- ceiling(max(b0, b1) / eqSolution$iStar)
          glrTables <- .updateTables(glrTables, b0=b0, b1=b1, p, N)
          estimate <- apply(glrTables$tauTable, 2, sum)
          names(estimate) <- c("alpha", "beta")
          if (verbose) {
            print("B0 going up")
            print(paste("Trying b0:", b0, "and b1:", b1))
            print(paste("New N", N))
            print("Estimates:")
            print(estimate)
          }
          resultList <- list(b1=c(resultList$b1, b1), b0=c(resultList$b0, b0),
                             estimate=rbind(resultList$estimate, estimate))
          if (estimate["alpha"] <= alpha) {
            break;
          }
        }
      }
      b0Limits[2] <- b0

    } else {

      b0Limits <- c(b0.1, b0)
      if (verbose) print("b0 values were banging...")
      ## change b1 after picking the b0 that goes below alpha
      if (estimate1["alpha"] <= alpha) {
        b1 <- b1.1
        b0 <- b0.1
        n2 <- n-1
      } # ow a and b remain the same..
      b1Limits <- c(b1, 0)
      if (estimate1["beta"] <= beta) {
        repeat {
          count <- count + 1
          if (count > maxIter) {
            ## We have exhausted our tries
            break;
          }
          b1 <-  b1 - stepSize
          N <- ceiling(max(b0, b1) / eqSolution$iStar)
          glrTables <- .updateTables(glrTables, b0=b0, b1=b1, p, N)
          estimate <- apply(glrTables$tauTable, 2, sum)
          names(estimate) <- c("alpha", "beta")
          if (verbose) {
            print("b1 going down")
            print(paste("Trying b0:", b0, "and b1:", b1))
            print(paste("New N", N))
            print("Estimates:")
            print(estimate)
          }

          resultList <- list(b1=c(resultList$b1, b1), b0=c(resultList$b0, b0),
                             estimate=rbind(resultList$estimate, estimate))
          if (estimate["beta"] >= beta) {
            break;
          }
        }
      } else {
        repeat {
          count <- count + 1
          if (count > maxIter) {
            ## We have exhausted our tries
            break;
          }
          b1 <-  b1 + stepSize
          N <- ceiling(max(b0, b1) / eqSolution$iStar)
          glrTables <- .updateTables(glrTables, b0=b0, b1=b1, p, N)
          estimate <- apply(glrTables$tauTable, 2, sum)
          names(estimate) <- c("alpha", "beta")
          if (verbose) {
            print("B0 going up")
            print(paste("Trying b0:", b0, "and b1:", b1))
            print(paste("New N", N))
            print("Estimates:")
            print(estimate)
          }
          resultList <- list(b1=c(resultList$b1, b1), b0=c(resultList$b0, b0),
                             estimate=rbind(resultList$estimate, estimate))
          if (estimate["beta"] <= beta) {
            break;
          }
        }
      }
      b1Limits[2] <- b1
    }
  } else {
      b1Limits <- resultList$b1[c(n-1, n)]
      b0Limits <- resultList$b0[c(n-1, n)]
  }

  ## Note down the number of iterations
  resultList$iterations <- count
  ##
  ## Now compute the values on a grid
  ##
  if (gridIt) {
    b1Limits <- sort(b1Limits)
    b0Limits <- sort(b0Limits)

    b1Vals <- seq(b1Limits[1], b1Limits[2], length.out=nGrid)
    b0Vals <- seq(b0Limits[1], b0Limits[2], length.out=nGrid)

    betaTable <- alphaTable <- matrix(0, nrow=nGrid, ncol=nGrid)
    rownames(alphaTable) <- rownames(betaTable) <- paste("b1=", round(b1Vals, 5),sep="")
    colnames(alphaTable) <- colnames(betaTable) <- paste("b0=", round(b0Vals, 5),sep="")

    for (i in 1:nGrid) {
      for(j in 1:nGrid) {
        if (verbose) {
          print(paste("Computing element:", i, j))
        }
        b1 <- b1Vals[i]
        b0 <- b0Vals[j]
        glrTables <- .updateTables(glrTables, b0=b0, b1=b1, p, N)
        estimate <- apply(glrTables$tauTable, 2, sum)
        names(estimate) <- c("alpha", "beta")
        alphaTable[i, j] <- estimate["alpha"]
        betaTable[i, j] <- estimate["beta"]
      }
    }
    resultList$alphaTable <- alphaTable
    resultList$betaTable <- betaTable
    resultList$b1Vals <- b1Vals
    resultList$b0Vals <- b0Vals
  }
  resultList$glrTables = glrTables
  return(resultList)
}

