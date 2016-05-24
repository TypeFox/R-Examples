
.computePStar <-
function(p, tol=1e-7) {
  ##
  ## This is the function used to compute the bounds
  ##
  eqn <- function(x, p) {
    p0 <- p[1]
    p1 <- p[2]
    if (x <= p0 || x >= p1) {
      ifelse (x >= p1, Inf, -Inf)
    } else {
      lhs <- x * log(x/p0) + (1-x) * log((1-x)/(1-p0))
      rhs <- x * log(x/p1) + (1-x) * log((1-x)/(1-p1))
      lhs - rhs
    }
  }
  p0 <- p[1]
  pStar <- uniroot(eqn, p, tol=tol, p=p)$root
  iStar <- pStar * log(pStar/p0) + (1-pStar) * log((1-pStar)/(1-p0))
  ##  print(paste("Pstar, iStar:", pStar, iStar))
  list(pStar=pStar, iStar=iStar)
}

## .expandListCapacity <-
## function(existingList, newLength) {
##   append(existingList, vector(mode="list", length=newLength - length(existingList)))
## }

.updateTables <-
function(glrTables, b0, b1, p, N) {

    g <- function(s, n, j) {
    pj <- p[j+1]
    ifelse(s==0,
           n * log(1/(1-pj)),
           ifelse(s == n,
                  n * log(1/pj),
                  { pHat <- s/n; qHat <- 1-pHat; s * log(pHat/pj) + (n-s) * log(qHat/(1-pj)) }
                  )
           )
  }    
  

  currentN <- length(glrTables$g0Table)
  if (N > currentN) {
    ##     g0Table <- .expandListCapacity(glrTables$g0Table, N)
    ##     g1Table <- .expandListCapacity(glrTables$g1Table, N)
    g0Table <- append(glrTables$g0Table, vector(mode = "list", length = N - currentN))
    g1Table <- append(glrTables$g1Table, vector(mode = "list", length = N - currentN))

    ## Now compute the new g values.
    for (n in (currentN+1):N) {
      ##print(n)
      g0Table[[n]] <- sapply(0:n, g, n=n, j=0)
      g1Table[[n]] <- sapply(0:n, g, n=n, j=1)      
    }
  } else {
    g0Table <- glrTables$g0Table
    g1Table <- glrTables$g1Table
  }

  p0Table <- p1Table <- vector(mode="list", length=N)
  
  p0 <- p[1]
  p1 <- p[2]
  
  for (n in 1:N) {
    p0Vec <- p1Vec <- numeric(n+1)
    names(p0Vec) <- names(p1Vec) <- as.character(0:n)
    
    g0Vec <- g0Table[[n]]
    g1Vec <- g1Table[[n]]
    
    for (s in 0:n) {
      g0 <- g0Vec[s+1]
      g1 <- g1Vec[s+1]
      
      if (n == 1) {
        p0Vec[s+1] <- ifelse(s==0, 1-p0, p0)
        p1Vec[s+1] <- ifelse(s==0, 1-p1, p1)
      } else {
        condition <- (g0 < b0) && (g1 < b1)
        p0Vec[s+1] <- ifelse(condition,
                             ifelse(s==0, 0, p0 * p0Table[[n-1]][s]) +
                             ifelse(s==n, 0, (1-p0) * p0Table[[n-1]][s + 1]),
                             0)
        p1Vec[s+1] <- ifelse(condition,
                             ifelse(s==0, 0, p1 * p1Table[[n-1]][s]) +
                             ifelse(s==n, 0, (1-p1) * p1Table[[n-1]][s + 1]),
                             0)
      }
    }

    p0Table[[n]] <- p0Vec
    p1Table[[n]] <- p1Vec
      
  }

  ## Compute \hat{alpha} and \hat{beta}, estimates of alpha and beta
  
  tauTable <- matrix(0.0, nrow=N, ncol=2)
  colnames(tauTable) <- c("b0", "b1")
  ## sbTable <- rep(0, N)
  ## slTable <- rep(N, N)
  for (n in 1:N) {
    startS <- ceiling(p0 * n)
    searchIndices <- (startS+1):(n+1)
    k <- which(g0Table[[n]][searchIndices] >= b0)
    if ((l <- length(k)) > 0) {
      sSmallest <- searchIndices[k[1]]
      ## sbTable[n] <- sSmallest
      ##       print(paste("p0: (n, k) = ", n, k))
      ##       print(paste("sSmallest:", sSmallest))
      ##       browser()
      tauTable[n, 1] <- ifelse(sSmallest == 1, 0, p0 * ifelse(n==1, 1, p0Table[[n-1]][sSmallest-1])) +
        ifelse(sSmallest == n+1, 0, (1-p0) * ifelse(n==1, 1, p0Table[[n-1]][sSmallest]))
    }
    
    endS <- floor(n * p1)
    k <- which(g1Table[[n]][1:(endS+1)] >= b1)
    
    if ((l <- length(k)) > 0) {
      sLargest <- k[l]
      ## slTable[n] <- sLargest
      ##       print(paste("p0: (n, k) = ", n, k))
      ##       print(paste("sLargest:", sLargest))      
      tauTable[n, 2] <- ifelse(sLargest == 1, 0, p1 * ifelse(n==1, 1, p1Table[[n-1]][sLargest-1])) +
        ifelse(sLargest == n+1, 0, (1-p1) * ifelse(n==1, 1, p1Table[[n-1]][sLargest]))
    }
  }
  
    list(g0Table=g0Table, g1Table=g1Table, tauTable=tauTable)
}

