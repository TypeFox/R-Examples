`modelChoose` <-
function(modelShow, criteria="BIC", tol=0.20) {
 x       <- modelShow
 maxCrit <- function(x) {
  if (any(is.na(x[(1:8), ]$SeT))) x <- x[-which(is.na(x[1:8, ]$SeT)),]
  maxCrit <- which(x[criteria] <= max(x[criteria], na.rm=TRUE)+tol & x[criteria] >= max(x[criteria], na.rm=TRUE)-tol)
  #if (any(is.na(x[maxCrit, ]$SeT))) maxCrit <- maxCrit[-which(is.na(x[maxCrit, ]$SeT))]
  #print(c(x[maxCrit, ]$ID, maxCrit, x[maxCrit, ]$LL, x[maxCrit, ]$SeT)) #
  minSeT  <- which(x[maxCrit,]$SeT == min(x[maxCrit,]$SeT, na.rm=TRUE))
  #if ( length(maxCrit) > 1) minSeT  <- which(x[maxCrit,]$SeT == min(x[maxCrit,]$SeT, na.rm=TRUE))
  #if (!length(maxCrit) > 1) minSeT  <- which(x[maxCrit,]$SeT == min(x[maxCrit,]$SeT, na.rm=TRUE))
  return(c(x[maxCrit[minSeT],]$MODEL))
  }
 return( by(x, list(ID=x$ID), maxCrit, simplify=TRUE) )
 }


