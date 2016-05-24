######################################################################
####functions for downstream analysis of mlDNA prediction results#####


.uniqueTissues <- function (x) {
    sampleNum <- ncol(x)
    if (is.null(colnames(x))) {
        tsMatrix <- diag(x = 1, nrow = sampleNum, ncol = sampleNum)
    }
    else {
        tsMatrix <- matrix(0, nrow = sampleNum, ncol = sampleNum)
        colnames(tsMatrix) <- colnames(x)
        uniTS <- c("")
        for (i in 1:sampleNum) {
            lastdot <- sapply(gregexpr("\\.", colnames(x)[i]), 
                tail, 1)
            if (lastdot < 0) {
                curTSName <- colnames(x)[i]
            }
            else {
                if (lastdot > 1) {
                  curTSName <- str_sub(colnames(x)[i], 1, lastdot - 
                    1)
                }
                else {
                  curTSName <- colnames(x)[i]
                }
            }
            curTSIndex <- which(uniTS == curTSName)
            if (length(curTSIndex) == 0) {
                if (i == 1) 
                  uniTS[1] <- curTSName
                else uniTS <- c(uniTS, curTSName)
                tsMatrix[length(uniTS), i] <- 1
            }
            else {
                tsMatrix[curTSIndex, i] <- 1
            }
        }
        tsMatrix <- tsMatrix[1:length(uniTS), ]
        rownames(tsMatrix) <- uniTS
    }
    tsMatrix
}



.getsgene <- function (x, Log = FALSE, Base = 2, AddOne = FALSE, tsThreshold = 0.95, MeanOrMax = "Mean", Fraction = TRUE) {
    if (AddOne) {
        x <- x + 1
    }else {
        subfun <- function(vec) {
            res <- FALSE
            if (sum(abs(vec)) == 0) 
                res <- TRUE
            res
        }
        zeroIndex <- which(apply(x, 1, subfun))
        if (length(zeroIndex) >= 1) 
            x <- x[-zeroIndex, ]
    }
    if (Log) {
        x <- log(x, Base)
    }
    onets <- function(vec, tsMatrix, MeanOrMax) {
        tscorematrix <- matrix(0, nrow = nrow(tsMatrix), ncol = 2)
        for (i in 1:nrow(tsMatrix)) {
            sampleIndex <- which(tsMatrix[i, ] > 0)
            meanvalue <- NULL
            if (MeanOrMax == "Mean") {
                meanvalue <- mean(vec[sampleIndex])
            }
            else {
                meanvalue <- max(vec[sampleIndex])
            }
            tscorematrix[i, 1] <- i
            if (meanvalue == 0) {
                meanvalue <- 1e-06
            }
            tscorematrix[i, 2] <- 1 - max(vec[-sampleIndex])/meanvalue
        }
        tmax = max(tscorematrix[, 2])
        tmaxidx = tscorematrix[which(tscorematrix[, 2] == tmax), 1][1]
        return(list(tmaxidx = tmaxidx, tmax = tmax, tscores = tscorematrix[,2] ))
    }


    tsMatrix <- uniqueTissues(x)
    tissueNum <- nrow(tsMatrix)
    tt <- apply(x, 1, onets, tsMatrix = tsMatrix, MeanOrMax = MeanOrMax)
    alltsmatrix <- matrix(0, nrow = nrow(x), ncol = tissueNum + 2 )
    rownames(alltsmatrix) <- rownames(x)
    colnames(alltsmatrix) <- c( rownames(tsMatrix), "CSMaxScore", "Condition" )
    for (i in 1:nrow(x)) {
        alltsmatrix[i, 1:tissueNum ] <- tt[[i]]$tscores
        alltsmatrix[i, tissueNum + 1] <- tt[[i]]$tmax
        alltsmatrix[i, tissueNum + 2] <- rownames(tsMatrix)[tt[[i]]$tmaxidx]
    }
    idx <- which(  as.numeric(alltsmatrix[, tissueNum + 1]) >= tsThreshold )
    tsgene <- NULL
    if( length(idx) >= 1 ) 
       tsgene <- rownames(x)[idx]

    if (Fraction) {
        #tsgene <- t(apply(tsgene, 1, function(vec) vec/sum(vec)))
    }
    return(list(alltscoremat = alltsmatrix, tsgene = tsgene, uniquets = tsMatrix))
}



ConditionSpecificGenes <- function( expmat, logtransformed = TRUE, base = 2, threshold = 0.75 ) {
  
    x <- expmat
    if( logtransformed  )
       x <- base^expmat
    res <- .getsgene(x = x, Log = FALSE, Base = base, AddOne = FALSE, tsThreshold = threshold, MeanOrMax = "Max", Fraction = FALSE)
    res <- list( CSGenes = res$tsgene, CSScoreMat = res$alltscoremat, uniqueCS = res$uniquets )
    res
}


ConvergenceDegree <- function( vec1, vec2 ) {
   
   if( length(vec1) == 0 | length(vec2) == 0 ) {
      stop("Error: no element in vec1 or vec2.\n")
   }
   length( intersect(vec1, vec2) )/min( length(vec1), length(vec2) )
}



AverageRankScore <- function( featureMat, selGenes ) {
  
  idx <- match( selGenes, rownames(featureMat) )
  m <- length(selGenes)
  n <- nrow(featureMat)
  
  res <- rep(0, ncol(featureMat))
  names(res) <- colnames(featureMat)
  for( i in 1:ncol(featureMat) ) {
    Rank <- rank( abs(featureMat[,i]) )
    RankSub <- Rank[idx]
    res[i] <- sum(RankSub)/(m*n)
  }#end for i
  
  res
}

