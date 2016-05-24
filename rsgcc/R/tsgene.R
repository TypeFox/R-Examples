##################################################################
##find tissue specific genes

if( !require(stringr) ) install.packages("stringr")
require(stringr)

uniqueTissues <- function(x) {
  
    sampleNum <- ncol(x)
    if( is.null(colnames(x)) ) {  #no annoted tissue sample, each sample belongs to different tissue by default
      tsMatrix <- diag(x = 1, nrow = sampleNum, ncol = sampleNum ) 
    }else {
      tsMatrix <- matrix(0, nrow = sampleNum, ncol = sampleNum )
      colnames(tsMatrix) <- colnames(x)
      uniTS <- c("")
      for( i in 1: sampleNum ) {
        lastdot <- sapply(gregexpr("\\.", colnames(x)[i]), tail, 1)
        if( lastdot < 0 )  { curTSName <- colnames(x)[i] 
        } else { 
          if( lastdot > 1 ) {curTSName <- str_sub( colnames(x)[i], 1, lastdot - 1) 
        }else {
          curTSName <- colnames(x)[i]
        }
      }
      
      curTSIndex <- which(uniTS == curTSName)
      if( length(curTSIndex) == 0 ) { #not record
        if( i == 1) uniTS[1] <- curTSName
        else        uniTS <- c(uniTS, curTSName)
        
        tsMatrix[length(uniTS), i] <- 1
      }else {
        tsMatrix[curTSIndex, i] <- 1
      } 
    }#end for i
    tsMatrix <- tsMatrix[1:length(uniTS),]
    rownames(tsMatrix) <- uniTS
  }#end else
    
    tsMatrix
  
}


##################################################################################
##ts gene: t measure, 1- max(nonTissue)/max(Tissue)
getsgene <- function(x, Log = FALSE, Base = 2, AddOne = FALSE, tsThreshold = 0.95, MeanOrMax = "Mean", Fraction = TRUE ) {
  
  ##remove all zeros
  if (AddOne) {
    x <- x + 1
  }else {
    subfun <- function( vec ) {
      res <- FALSE
      if( sum(abs(vec)) == 0 ) 
        res <- TRUE
      res
    }
    zeroIndex <- which(apply(x, 1, subfun ) )
    if( length(zeroIndex) >= 1 )
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
      if( MeanOrMax == "Mean" ) {
         meanvalue <- mean(vec[sampleIndex])  ##consider the max values in one stress
       }else {
         meanvalue <- max(vec[sampleIndex])  ##consider the max values in one stress
       }
      tscorematrix[i, 1] <- i
      if (meanvalue == 0) {
        meanvalue <- 1.0E-6
      }
      tscorematrix[i, 2] <- 1 - max(vec[-sampleIndex])/meanvalue
    }
    tmax = max(tscorematrix[, 2])
    tmaxidx = tscorematrix[which(tscorematrix[, 2] == tmax), 1][1]
    return(list(tmaxidx = tmaxidx, tmax = tmax, tscorematrix = tscorematrix))
  }
  


  tsMatrix <- uniqueTissues(x)
  tt <- apply(x, 1, onets, tsMatrix = tsMatrix, MeanOrMax = MeanOrMax)
  
  alltsmatrix <- matrix( 0, nrow = nrow(x), ncol = nrow(tsMatrix) )
  rownames(alltsmatrix) <- rownames(x)
  colnames(alltsmatrix) <- rownames(tsMatrix)
  tscorematrix <- matrix(0, nrow = nrow(x), ncol = 3)
  rownames(tscorematrix) <- rownames(x)
  colnames(tscorematrix) <- c("GeneIndex", "tsmaxscore", "tsmaxidx")
  for (i in 1:nrow(x)) {
    tscorematrix[i, 1] <- i
    tscorematrix[i, 2] <- tt[[i]]$tmax
    tscorematrix[i, 3] <- tt[[i]]$tmaxidx 
    alltsmatrix[i,] <- tt[[i]]$tscorematrix[,2]
  }
  
  
  tscore <- tscorematrix[which(tscorematrix[, 2] >= tsThreshold), ]
  tsgene <- x[tscore[, 1], ]
  
  if (Fraction) {
    tsgene <- t(apply(tsgene, 1, function(vec) vec/sum(vec)))
  }
  
  return(list(alltscoremat = alltsmatrix,  tscore = tscorematrix, tsgene = tsgene, uniquets = tsMatrix))

}
 
