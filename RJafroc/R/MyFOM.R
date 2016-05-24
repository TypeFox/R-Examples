MyFOM <- function(nl, ll, lesionNum, lesionID, lesionWeight, maxNL, fom) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  K <- length(nl[, 1])
  K2 <- length(ll[, 1])
  K1 <- K - K2
  
  FOM <- NA
  if (fom == "Wilcoxon") {
    truth <- c(rep(0, K1), rep(1, K2))
    ratings <- c(nl[1:K1], ll)
    FOM <- CalculateTrapezoidalArea(ratings, truth)
  } else if (fom == "HrAuc") {
    truth <- c(rep(0, K1), rep(1, K2))
    ratings <- array(dim = c(K1 + K2))
    
    for (k in 1:K1) {
      ratings[k] <- max(nl[k, ])
    }
    
    
    for (k in (K1 + 1):(K1 + K2)) {
      ratings[k] <- max(c(nl[k, ], ll[k - K1, ]))
    }
    
    FOM <- CalculateTrapezoidalArea(ratings, truth)
  } else if (fom == "HrSe") {
    tpCount <- 0
    for (k in 1:K2) {
      tempMax <- UNINITIALIZED
      for (el in 1:lesionNum[k]) tempMax <- max(tempMax, ll[k, el])
      
      for (el in 1:maxNL) tempMax <- max(tempMax, nl[k + K1, el])
      
      if (tempMax > UNINITIALIZED) 
        tpCount <- tpCount + 1
    }
    FOM <- tpCount/K2
  } else if (fom == "HrSp") {
    fpCount <- 0
    for (k in 1:K1) {
      tempMax <- UNINITIALIZED
      for (el in 1:maxNL) tempMax <- max(tempMax, nl[k, el])
      if (tempMax > UNINITIALIZED) 
        fpCount <- fpCount + 1
    }
    FOM <- 1 - fpCount/K1
  } else if (fom == "SongA1") {
    x <- array(UNINITIALIZED, dim = c(K))
    y1 <- array(UNINITIALIZED, dim = c(K2))
    y2 <- array(UNINITIALIZED, dim = c(K2))
    
    kInd <- which(apply(nl != UNINITIALIZED, 1, any))
    for (k in kInd){      
      x[ k ] <- mean(nl[ k, ( nl[ k, ] != UNINITIALIZED )])
    }
    y1 <- x[(K1 + 1):K]
    x <- x[1:K1]
    px0 <- sum(x == UNINITIALIZED)    
    
    kInd <- which(apply(ll != UNINITIALIZED, 1, any))
    for (k in kInd){      
      y2[ k ] <- mean(ll[ k, ( ll[ k, ] != UNINITIALIZED )])
    }
    y <- apply(cbind(y1, y2), 1, max)
    py0 <- sum(y == UNINITIALIZED)
    
    FOM <- 0
    x <- x[x != UNINITIALIZED]
    y <- y[y != UNINITIALIZED]
    for ( k1 in 1 : length(x)){
      FOM <- FOM + sum(x[k1] < y) + 0.5 * sum(x[k1] == y)
    }
    FOM <- FOM / (K1 * K2) + (px0 / K1) * ( 1 - 0.5 * py0 / K2)
  } else if (fom == "SongA2") {
    n0 <- apply(nl[1:K1, ] != UNINITIALIZED, 1, sum)
    n1 <- apply(cbind(nl[(K1 + 1):K, ], ll) != UNINITIALIZED, 1, sum)   
    k1Ind <- which(n0 > 0)
    k2Ind <- which(n1 > 0)
    px0 <- K1 - length(k1Ind)
    py0 <- K2 - length(k2Ind)
    v1 <- 0.5
    tw <- 0
    for (k1 in k1Ind) {
      for (k2 in k2Ind) {          
        v2 <- 0
        x <- nl[k1, ][nl[k1, ] != UNINITIALIZED]
        y1 <- nl[K1 + k2, ][nl[K1 + k2, ] != UNINITIALIZED]
        y2 <- ll[k2, ][ll[k2, ] != UNINITIALIZED]
        for (l in 1:length(x)){
          v2 <- v2 + sum(x[l] < y1) + 0.5 * sum(x[l] == y1) + sum(x[l] < y2) + 0.5 * sum(x[l] == y2)
        }          
        v2 <- v2/(n0[k1] * n1[k2])
        tw <- tw + (v1 < v2) + 0.5 * (v1 == v2)
      }
    }
    FOM <- tw/(K1 * K2) + px0/K1 * (1 - 0.5 * py0/K2)
  } else if (fom == "JAFROC1") {
    numLesTotal <- sum(lesionNum)
    les <- as.vector(ll[lesionID != UNINITIALIZED])
    fp <- array(dim = c(K))
    for (k in 1:K) {
      fp[k] <- max(nl[k, ])
    }
    ratings <- c(fp, les)
    truth <- c(rep(0, K), rep(1, numLesTotal))
    FOM <- CalculateTrapezoidalArea(ratings, truth)
  } else if (fom == "JAFROC") {
    numLesTotal <- sum(lesionNum)
    les <- as.vector(ll[lesionID != UNINITIALIZED])
    fp <- array(dim = c(K1))
    for (k in 1:K1) {
      fp[k] <- max(nl[k, ])
    }
    ratings <- c(fp, les)
    truth <- c(rep(0, K1), rep(1, length(les)))
    FOM <- CalculateTrapezoidalArea(ratings, truth)
  } else if (fom == "wJAFROC1") {
    numLesTotal <- sum(lesionNum)
    les <- as.vector(ll[lesionID != UNINITIALIZED])
    fp <- array(dim = c(K))
    for (k in 1:K) {
      fp[k] <- max(nl[k, ])
    }
    ratings <- c(fp, les)
    truth <- c(rep(0, K), rep(1, numLesTotal))
    weights <- lesionWeight[lesionWeight != UNINITIALIZED]
    FOM <- CalculateTrapezoidalAreaWeighted(ratings, truth, weights, K2)
  } else if (fom == "wJAFROC") {
    numLesTotal <- sum(lesionNum)
    les <- as.vector(ll[lesionID != UNINITIALIZED])
    fp <- array(dim = c(K1))
    for (k in 1:K1) {
      fp[k] <- max(nl[k, ])
    }
    ratings <- c(fp, les)
    truth <- c(rep(0, K1), rep(1, numLesTotal))
    weights <- lesionWeight[lesionWeight != UNINITIALIZED]
    FOM <- CalculateTrapezoidalAreaWeighted(ratings, truth, weights, K2)
  } else if (fom == "MaxLLF") {
    numMarksTotal <- sum(ll != UNINITIALIZED)
    numLesTotal <- sum(lesionNum)
    FOM <- numMarksTotal/numLesTotal
  } else if (fom == "MaxNLF") {
    numMarksTotal <- sum(nl[1:K1, ] != UNINITIALIZED)
    FOM <- numMarksTotal/K1
  } else if (fom == "MaxNLFAllCases") {
    numMarksTotal <- sum(nl != UNINITIALIZED)
    FOM <- numMarksTotal/K
  } else if (fom == "ExpTrnsfmSp") {
    numMarksTotal <- sum(nl[1:K1, ] != UNINITIALIZED)
    FOM <- exp(-numMarksTotal/K1)
  } else if (fom == "ROI") {
    nn <- 0
    ns <- sum(lesionNum)
    tw <- 0
    for (k in 1:K) {
      for (el in 1:maxNL) {
        if (nl[k, el] == UNINITIALIZED) 
          next
        nn <- nn + 1
        tw <- tw + sum(nl[k, el] < as.vector(ll)) + 0.5 * sum(nl[k, el] == as.vector(ll))
      }
    }
    FOM <- tw/(nn * ns)
  } else {
    errMsg <- paste0(fom, " is not an available figure of merit.")
    stop(errMsg)
  }
  return(FOM)
} 
