##summarize detection histories and count data
countDist <- function(object, plot.freq = TRUE, plot.distance = TRUE, ...){
  UseMethod("countDist", object)
}



countDist.default <- function(object, plot.freq = TRUE, plot.distance = TRUE, ...){
  stop("\nFunction not yet defined for this object class\n")
}




##for unmarkedFrameDS
countDist.unmarkedFrameDS <- function(object, plot.freq = TRUE, plot.distance = TRUE, ...) {

  ##extract data
  yMat <- object@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  nvisits <- 1

  ##distance classes
  dist.classes <- object@dist.breaks

  ##number of distance classes
  n.dist.classes <- length(dist.classes) - 1
    
  ##units
  unitsIn <- object@unitsIn
  
  ##create string of names
  dist.names <- rep(NA, n.dist.classes)
  for(i in 1:n.dist.classes){
    dist.names[i] <- paste(dist.classes[i], "-", dist.classes[i+1], sep = "")
  }
  
  ##visits per season
  n.visits.season <- nvisits/n.seasons

  ##collapse yMat into a single vector
  yVec <- as.vector(yMat)

  ##determine size of plot window
  ##when both types are requested
  if(plot.freq && plot.distance) {
    par(mfrow = c(1, 2))
  }
  
  ##summarize counts
  if(plot.freq) {
    
    ##create histogram
    barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
         main = "Distribution of raw counts", cex.lab = 1.2)
  }

  ##summarize counts per distance
  dist.sums.full <- colSums(yMat)
  names(dist.sums.full) <- dist.names
  dist.sums.seasons <- list(dist.sums.full)
  
  if(plot.distance) {
    
    ##create histogram
    barplot(dist.sums.full, ylab = "Frequency",
            xlab = paste("Distance class (", unitsIn, ")", sep = ""),
            main = "Distribution of distance data", cex.lab = 1.2)
  }

  ##raw counts
  count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
  count.table.seasons <- list(count.table.full)


  ##for each season, determine frequencies
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  ySeason <- yMat
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(ySeason, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
  ##number of sites sampled
  out.freqs[1, 1] <- sum(!is.na(sum.rows))
  ##number of sites with at least 1 detection
  out.freqs[1, 2] <- sum(det.sum)

  ##create a matrix with proportion of sites with colonizations
  ##and extinctions based on raw data
  out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
  colnames(out.props) <- "naive.occ"
  rownames(out.props) <- rownames(out.freqs)
  out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]

  out.count <- list("count.table.full" = count.table.full,
                    "count.table.seasons" = count.table.seasons,
                    "dist.sums.full" = dist.sums.full,
                    "dist.sums.seasons" = dist.sums.seasons,
                    "dist.names" = dist.names,
                    "n.dist.classes" = n.dist.classes,
                    "out.freqs" = out.freqs,
                    "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season)
  class(out.count) <- "countDist"
  return(out.count)
}



##for unmarkedFitDS
countDist.unmarkedFitDS <- function(object, plot.freq = TRUE, plot.distance = TRUE, ...) {

  ##extract data
  yMat <- object@data@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  nvisits <- 1

  ##distance classes
  dist.classes <- object@data@dist.breaks

  ##number of distance classes
  n.dist.classes <- length(dist.classes) - 1

  ##units
  unitsIn <- object@data@unitsIn
  
  ##create string of names
  dist.names <- rep(NA, n.dist.classes)
  for(i in 1:n.dist.classes){
    dist.names[i] <- paste(dist.classes[i], "-", dist.classes[i+1], sep = "")
  }
  
  ##visits per season
  n.visits.season <- nvisits/n.seasons

  ##collapse yMat into a single vector
  yVec <- as.vector(yMat)

  ##determine size of plot window
  ##when both types are requested
  if(plot.freq && plot.distance) {
    par(mfrow = c(1, 2))
  }
  
  ##summarize counts
  if(plot.freq) {
    
    ##create histogram
    barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
         main = "Distribution of raw counts", cex.lab = 1.2)
  }

  ##summarize counts per distance
  dist.sums.full <- colSums(yMat)
  names(dist.sums.full) <- dist.names
  dist.sums.seasons <- list(dist.sums.full)
  
  if(plot.distance) {
    
    ##create histogram
    barplot(dist.sums.full, ylab = "Frequency",
            xlab = paste("Distance class (", unitsIn, ")", sep = ""),
            main = "Distribution of distance data", cex.lab = 1.2)
  }

  ##raw counts
  count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
  count.table.seasons <- list(count.table.full)


  ##for each season, determine frequencies
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  ySeason <- yMat
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(ySeason, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
  ##number of sites sampled
  out.freqs[1, 1] <- sum(!is.na(sum.rows))
  ##number of sites with at least 1 detection
  out.freqs[1, 2] <- sum(det.sum)

  ##create a matrix with proportion of sites with colonizations
  ##and extinctions based on raw data
  out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
  colnames(out.props) <- "naive.occ"
  rownames(out.props) <- rownames(out.freqs)
  out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]

  out.count <- list("count.table.full" = count.table.full,
                    "count.table.seasons" = count.table.seasons,
                    "dist.sums.full" = dist.sums.full,
                    "dist.sums.seasons" = dist.sums.seasons,
                    "dist.names" = dist.names,
                    "n.dist.classes" = n.dist.classes,
                    "out.freqs" = out.freqs,
                    "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season)
  class(out.count) <- "countDist"
  return(out.count)
}



##for unmarkedFrameGDS
countDist.unmarkedFrameGDS <- function(object, plot.freq = TRUE, plot.distance = TRUE, ...) {

  ##extract data
  yMat <- object@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  ##for GDS - several visits in single season
  nvisits <- object@numPrimary

  ##distance classes
  dist.classes <- object@dist.breaks

  ##number of distance classes
  n.dist.classes <- length(dist.classes) - 1

  ##units
  unitsIn <- object@unitsIn
  
  ##create string of names
  dist.names <- rep(NA, n.dist.classes)
  for(i in 1:n.dist.classes){
    dist.names[i] <- paste(dist.classes[i], "-", dist.classes[i+1], sep = "")
  }
  
  ##visits per season
  n.visits.season <- nvisits/n.seasons

  ##collapse yMat into a single vector
  yVec <- as.vector(yMat)

  ##determine size of plot window
  ##when both types are requested
  if(plot.freq && plot.distance) {
    par(mfrow = c(1, 2))
  }
  
  ##summarize counts
  if(plot.freq) {
    
    ##create histogram
    barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
         main = "Distribution of raw counts", cex.lab = 1.2)
  }

  ##summarize counts per distance
  dist.sums.full <- rep(NA, n.dist.classes)
  ##create matrix to hold indices of visits x dist.classes
  mat.dist <- matrix(1:(nvisits*n.dist.classes),
                     nrow = nvisits,
                     ncol = n.dist.classes)
  for(j in 1:n.dist.classes) {
    dist.sums.full[j] <- sum(colSums(yMat[, mat.dist[j, ]]))
  }
  names(dist.sums.full) <- dist.names
  dist.sums.seasons <- list(dist.sums.full)
  
  if(plot.distance) {
    
    ##create histogram
    barplot(dist.sums.full, ylab = "Frequency",
            xlab = paste("Distance class (", unitsIn, ")", sep = ""),
            main = "Distribution of distance data", cex.lab = 1.2)
  }

  ##raw counts
  count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
  count.table.seasons <- list(count.table.full)


  ##for each season, determine frequencies
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  ySeason <- yMat
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(ySeason, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
  ##number of sites sampled
  out.freqs[1, 1] <- sum(!is.na(sum.rows))
  ##number of sites with at least 1 detection
  out.freqs[1, 2] <- sum(det.sum)

  ##create a matrix with proportion of sites with colonizations
  ##and extinctions based on raw data
  out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
  colnames(out.props) <- "naive.occ"
  rownames(out.props) <- rownames(out.freqs)
  out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]

  out.count <- list("count.table.full" = count.table.full,
                    "count.table.seasons" = count.table.seasons,
                    "dist.sums.full" = dist.sums.full,
                    "dist.sums.seasons" = dist.sums.seasons,
                    "dist.names" = dist.names,
                    "n.dist.classes" = n.dist.classes,
                    "out.freqs" = out.freqs,
                    "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season)
  class(out.count) <- "countDist"
  return(out.count)
}



##for unmarkedFitGDS
countDist.unmarkedFitGDS <- function(object, plot.freq = TRUE, plot.distance = TRUE, ...) {

  ##extract data
  yMat <- object@data@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  ##for GDS - several visits in single season
  nvisits <- object@data@numPrimary

  ##distance classes
  dist.classes <- object@data@dist.breaks

  ##number of distance classes
  n.dist.classes <- length(dist.classes) - 1

  ##units
  unitsIn <- object@data@unitsIn
  
  ##create string of names
  dist.names <- rep(NA, n.dist.classes)
  for(i in 1:n.dist.classes){
    dist.names[i] <- paste(dist.classes[i], "-", dist.classes[i+1], sep = "")
  }
  
  ##visits per season
  n.visits.season <- nvisits/n.seasons

  ##collapse yMat into a single vector
  yVec <- as.vector(yMat)

  ##determine size of plot window
  ##when both types are requested
  if(plot.freq && plot.distance) {
    par(mfrow = c(1, 2))
  }
  
  ##summarize counts
  if(plot.freq) {
    
    ##create histogram
    barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
         main = "Distribution of raw counts", cex.lab = 1.2)
  }

  ##summarize counts per distance
  dist.sums.full <- rep(NA, n.dist.classes)
  ##create matrix to hold indices of visits x dist.classes
  mat.dist <- matrix(1:(nvisits*n.dist.classes),
                     nrow = nvisits,
                     ncol = n.dist.classes)
  for(j in 1:n.dist.classes) {
    dist.sums.full[j] <- sum(colSums(yMat[, mat.dist[j, ]]))
  }
  names(dist.sums.full) <- dist.names
  dist.sums.seasons <- list(dist.sums.full)
  
  if(plot.distance) {
    
    ##create histogram
    barplot(dist.sums.full, ylab = "Frequency",
            xlab = paste("Distance class (", unitsIn, ")", sep = ""),
            main = "Distribution of distance data", cex.lab = 1.2)
  }

  ##raw counts
  count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
  count.table.seasons <- list(count.table.full)


  ##for each season, determine frequencies
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  ySeason <- yMat
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(ySeason, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
  ##number of sites sampled
  out.freqs[1, 1] <- sum(!is.na(sum.rows))
  ##number of sites with at least 1 detection
  out.freqs[1, 2] <- sum(det.sum)

  ##create a matrix with proportion of sites with colonizations
  ##and extinctions based on raw data
  out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
  colnames(out.props) <- "naive.occ"
  rownames(out.props) <- rownames(out.freqs)
  out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]

  out.count <- list("count.table.full" = count.table.full,
                    "count.table.seasons" = count.table.seasons,
                    "dist.sums.full" = dist.sums.full,
                    "dist.sums.seasons" = dist.sums.seasons,
                    "dist.names" = dist.names,
                    "n.dist.classes" = n.dist.classes,
                    "out.freqs" = out.freqs,
                    "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season)
  class(out.count) <- "countDist"
  return(out.count)
}



##print method
print.countDist <- function(x, digits = 2, ...) {
  #if(x$n.seasons == 1) {
    cat("\nSummary of counts:\n")
    count.mat <- matrix(x$count.table.full, nrow = 1)
    colnames(count.mat) <- names(x$count.table.full)
    rownames(count.mat) <- "Frequency"
    print(count.mat)
    
    cat("\nSummary of distance data:\n")
    out.mat <- matrix(x$dist.sums.full, nrow = 1)
    colnames(out.mat) <- names(x$dist.sums.full)
    rownames(out.mat) <- "Frequency"
    print(out.mat)

    cat("\nProportion of sites with at least one detection:\n", round(x$out.props[, "naive.occ"], digits), "\n\n")
    
    cat("Frequencies of sites with detections:\n")
    ##add matrix of frequencies
    print(x$out.freqs)
#  }

    ##for data across several seasons, display in a matrix
    ##det.dist <- matrix(unlist(x$dist.sums.seasons),
    ##                   nrow = x$n.seasons)
    ##colnames(det.dist) <- names(x$dist.sums.full)
    ##rownames(det.dist) <- paste("Season-", 1:x$n.seasons, sep = "")
  }    

