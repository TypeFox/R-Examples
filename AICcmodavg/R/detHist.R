##summarize detection histories and count data
detHist <- function(object, ...){
  UseMethod("detHist", object)
}



detHist.default <- function(object, ...){
  stop("\nFunction not yet defined for this object class\n")
}



##for unmarkedFrameOccu (same as data format for occuRN)
detHist.unmarkedFrameOccu <- function(object, ...) {

  ##extract data
  yMat <- object@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  nvisits <- ncol(yMat)
  ##visits per season
  n.visits.season <- nvisits/n.seasons
  
  ##summarize detection histories
  hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
  hist.table.full <- table(hist.full, deparse.level = 0)

  ##for each season, determine frequencies
  hist.table.seasons <- vector(mode = "list", length = n.seasons)
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  hist.table.seasons[[1]]$hist.table <- hist.table.full
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = yMat, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(yMat, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(yMat)) == ncol(yMat)
    
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

  out.det <- list("hist.table.full" = hist.table.full,
                  "hist.table.seasons" = hist.table.seasons,
                  "out.freqs" = out.freqs, "out.props" = out.props,
                  "n.seasons" = n.seasons,
                  "n.visits.season" = n.visits.season)
  class(out.det) <- "detHist"
  return(out.det)
}



##for occu
detHist.unmarkedFitOccu <- function(object, ...) {

  ##extract data
  yMat <- object@data@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  nvisits <- ncol(yMat)
  ##visits per season
  n.visits.season <- nvisits/n.seasons
  
  ##summarize detection histories
  hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
  hist.table.full <- table(hist.full, deparse.level = 0)

  ##for each season, determine frequencies
  hist.table.seasons <- vector(mode = "list", length = n.seasons)
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  hist.table.seasons[[1]]$hist.table <- hist.table.full
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = yMat, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(yMat, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(yMat)) == ncol(yMat)
    
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

  out.det <- list("hist.table.full" = hist.table.full,
                  "hist.table.seasons" = hist.table.seasons,
                  "out.freqs" = out.freqs, "out.props" = out.props,
                  "n.seasons" = n.seasons,
                  "n.visits.season" = n.visits.season)
  class(out.det) <- "detHist"
  return(out.det)
}



##for unmarkedFrameOccuFP
detHist.unmarkedFrameOccuFP <- function(object, ...) {

  ##extract data
  yMat <- object@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  nvisits <- ncol(yMat)
  ##visits per season
  n.visits.season <- nvisits/n.seasons
  
  ##summarize detection histories
  hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
  hist.table.full <- table(hist.full, deparse.level = 0)

  ##for each season, determine frequencies
  hist.table.seasons <- vector(mode = "list", length = n.seasons)
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  hist.table.seasons[[1]]$hist.table <- hist.table.full
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = yMat, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(yMat, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(yMat)) == ncol(yMat)
    
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

  out.det <- list("hist.table.full" = hist.table.full,
                  "hist.table.seasons" = hist.table.seasons,
                  "out.freqs" = out.freqs, "out.props" = out.props,
                  "n.seasons" = n.seasons,
                  "n.visits.season" = n.visits.season)
  class(out.det) <- "detHist"
  return(out.det)
}



##for occuFP
detHist.unmarkedFitOccuFP <- function(object, ...) {

  ##extract data
  yMat <- object@data@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  nvisits <- ncol(yMat)
  ##visits per season
  n.visits.season <- nvisits/n.seasons
  
  ##summarize detection histories
  hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
  hist.table.full <- table(hist.full, deparse.level = 0)

  ##for each season, determine frequencies
  hist.table.seasons <- vector(mode = "list", length = n.seasons)
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  hist.table.seasons[[1]]$hist.table <- hist.table.full
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = yMat, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(yMat, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(yMat)) == ncol(yMat)
    
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

  out.det <- list("hist.table.full" = hist.table.full,
                  "hist.table.seasons" = hist.table.seasons,
                  "out.freqs" = out.freqs, "out.props" = out.props,
                  "n.seasons" = n.seasons,
                  "n.visits.season" = n.visits.season)
  class(out.det) <- "detHist"
  return(out.det)
}



##for occuRN
detHist.unmarkedFitOccuRN <- function(object, ...) {

  ##extract data
  yMat <- object@data@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  nvisits <- ncol(yMat)
  ##visits per season
  n.visits.season <- nvisits/n.seasons
  
  ##summarize detection histories
  hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
  hist.table.full <- table(hist.full, deparse.level = 0)

  ##for each season, determine frequencies
  hist.table.seasons <- vector(mode = "list", length = n.seasons)
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  hist.table.seasons[[1]]$hist.table <- hist.table.full
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = yMat, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(yMat, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(yMat)) == ncol(yMat)
    
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

  out.det <- list("hist.table.full" = hist.table.full,
                  "hist.table.seasons" = hist.table.seasons,
                  "out.freqs" = out.freqs, "out.props" = out.props,
                  "n.seasons" = n.seasons,
                  "n.visits.season" = n.visits.season)
  class(out.det) <- "detHist"
  return(out.det)
}



##for unmarkedMultFrame
detHist.unmarkedMultFrame <- function(object, ...) {

  ##extract data
  yMat <- object@y
  nsites <- nrow(yMat)
  n.seasons <- object@numPrimary
  nvisits <- ncol(yMat)
  ##visits per season
  n.visits.season <- nvisits/n.seasons

  ##summarize detection histories
  hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
  hist.table.full <- table(hist.full, deparse.level = 0)

  ##for each season, determine frequencies
  out.seasons <- vector(mode = "list", length = n.seasons)
  hist.table.seasons <- vector(mode = "list", length = n.seasons)
  out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected", "colonized",
                           "extinct", "static", "common")
  rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")

  ##sequence of visits
  vis.seq <- seq(from = 1, to = nvisits, by = n.visits.season)
  for(i in 1:n.seasons) {
    col.start <- vis.seq[i]
    col.end <- col.start + (n.visits.season - 1)
    ySeason <- yMat[, col.start:col.end]
    ##summarize detection histories
    det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
    hist.table.seasons[[i]]$hist.table <- table(det.hist, deparse.level = 0)

    ##determine proportion of sites with at least 1 detection
    det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

    ##check sites with observed detections and deal with NA's
    sum.rows <- rowSums(ySeason, na.rm = TRUE)
    is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
    ##number of sites sampled
    out.freqs[i, 1] <- sum(!is.na(sum.rows))
    out.freqs[i, 2] <- sum(det.sum)


    #sites without detections
    none <- which(sum.rows == 0)
    #sites with at least one detection
    some <- which(sum.rows != 0) 
    out.seasons[[i]] <- list("none" = none, "some" = some)
  }
  
  ##populate out.freqs with freqs of extinctions and colonizations
  for(j in 2:n.seasons) {
    none1 <- out.seasons[[j-1]]$none
    some1 <- out.seasons[[j-1]]$some
    none2 <- out.seasons[[j]]$none
    some2 <- out.seasons[[j]]$some
    ##colonizations
    out.freqs[j, 3] <- sum(duplicated(c(some2, none1)))
    ##extinctions
    out.freqs[j, 4] <- sum(duplicated(c(some1, none2)))
    ##no change
    out.freqs[j, 5] <- sum(duplicated(c(some1, some2))) + sum(duplicated(c(none1, none2)))
    ##sites both sampled in t and t-1
    year1 <- c(none1, some1)
    year2 <- c(none2, some2)
    out.freqs[j, 6] <- sum(duplicated(c(year1, year2)))
  }

  ##create a matrix with proportion of sites with colonizations
  ##and extinctions based on raw data
  out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 4)
  colnames(out.props) <- c("naive.occ", "naive.colonization", "naive.extinction", "naive.static")
  rownames(out.props) <- rownames(out.freqs)
  out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]
  out.props[, 2] <- out.freqs[, 3]/out.freqs[, 6]
  out.props[, 3] <- out.freqs[, 4]/out.freqs[, 6]
  out.props[, 4] <- out.freqs[, 5]/out.freqs[, 6]
  
  out.det <- list("hist.table.full" = hist.table.full,
                  "hist.table.seasons" = hist.table.seasons,
                  "out.freqs" = out.freqs, "out.props" = out.props,
                  "n.seasons" = n.seasons,
                  "n.visits.season" = n.visits.season)
  class(out.det) <- "detHist"
  return(out.det)
}



##for colext
detHist.unmarkedFitColExt <- function(object, ...) {

  ##extract data
  yMat <- object@data@y
  nsites <- nrow(yMat)
  n.seasons <- object@data@numPrimary
  nvisits <- ncol(yMat)
  ##visits per season
  n.visits.season <- nvisits/n.seasons

  ##summarize detection histories
  hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
  hist.table.full <- table(hist.full, deparse.level = 0)

  ##for each season, determine frequencies
  out.seasons <- vector(mode = "list", length = n.seasons)
  hist.table.seasons <- vector(mode = "list", length = n.seasons)
  out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected", "colonized",
                           "extinct", "static", "common")
  rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")

  ##sequence of visits
  vis.seq <- seq(from = 1, to = nvisits, by = n.visits.season)
  for(i in 1:n.seasons) {
    col.start <- vis.seq[i]
    col.end <- col.start + (n.visits.season - 1)
    ySeason <- yMat[, col.start:col.end]
    ##summarize detection histories
    det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
    hist.table.seasons[[i]]$hist.table <- table(det.hist, deparse.level = 0)

    ##determine proportion of sites with at least 1 detection
    det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

    ##check sites with observed detections and deal with NA's
    sum.rows <- rowSums(ySeason, na.rm = TRUE)
    is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
    ##number of sites sampled
    out.freqs[i, 1] <- sum(!is.na(sum.rows))
    out.freqs[i, 2] <- sum(det.sum)


    #sites without detections
    none <- which(sum.rows == 0)
    #sites with at least one detection
    some <- which(sum.rows != 0) 
    out.seasons[[i]] <- list("none" = none, "some" = some)
  }
  
  ##populate out.freqs with freqs of extinctions and colonizations
  for(j in 2:n.seasons) {
    none1 <- out.seasons[[j-1]]$none
    some1 <- out.seasons[[j-1]]$some
    none2 <- out.seasons[[j]]$none
    some2 <- out.seasons[[j]]$some
    ##colonizations
    out.freqs[j, 3] <- sum(duplicated(c(some2, none1)))
    ##extinctions
    out.freqs[j, 4] <- sum(duplicated(c(some1, none2)))
    ##no change
    out.freqs[j, 5] <- sum(duplicated(c(some1, some2))) + sum(duplicated(c(none1, none2)))
    ##sites both sampled in t and t-1
    year1 <- c(none1, some1)
    year2 <- c(none2, some2)
    out.freqs[j, 6] <- sum(duplicated(c(year1, year2)))
  }

  ##create a matrix with proportion of sites with colonizations
  ##and extinctions based on raw data
  out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 4)
  colnames(out.props) <- c("naive.occ", "naive.colonization", "naive.extinction", "naive.static")
  rownames(out.props) <- rownames(out.freqs)
  out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]
  out.props[, 2] <- out.freqs[, 3]/out.freqs[, 6]
  out.props[, 3] <- out.freqs[, 4]/out.freqs[, 6]
  out.props[, 4] <- out.freqs[, 5]/out.freqs[, 6]
  
  out.det <- list("hist.table.full" = hist.table.full,
                  "hist.table.seasons" = hist.table.seasons,
                  "out.freqs" = out.freqs, "out.props" = out.props,
                  "n.seasons" = n.seasons,
                  "n.visits.season" = n.visits.season)
  class(out.det) <- "detHist"
  return(out.det)
}



##print method
print.detHist <- function(x, digits = 2, ...) {
  if(identical(x$n.seasons, 1)) {
    cat("\nSummary of detection histories: \n")
    out.mat <- matrix(x$hist.table.full, nrow = 1)
    colnames(out.mat) <- names(x$hist.table.full)
    rownames(out.mat) <- "Frequency"
    print(out.mat)
    cat("\nProportion of sites with at least one detection:\n", round(x$out.props[, "naive.occ"], digits), "\n\n")

    cat("Frequencies of sites with detections:\n")
    ##add matrix of frequencies
    print(x$out.freqs)

  } else {
    cat("\nSummary of detection histories (", x$n.seasons, " seasons combined): \n", sep ="")
    ##determine number of characters
    num.chars <- nchar(paste(names(x$hist.table.full), collapse = ""))
    if(num.chars >= 80) {
      cat("\nNote:  Detection histories exceed 80 characters and are not displayed\n")
    } else {
        out.mat <- matrix(x$hist.table.full, nrow = 1)
        colnames(out.mat) <- names(x$hist.table.full)
        rownames(out.mat) <- "Frequency"
        print(out.mat)
      }
    
      cat("\nSeason-specific detection histories: \n")
      for(i in 1:x$n.seasons) {
        cat("Season", i, "\n")
        temp.tab <- x$hist.table.seasons[[i]]$hist.table
        out.mat <- matrix(temp.tab, nrow = 1)
        colnames(out.mat) <- names(temp.tab)
        rownames(out.mat) <- "Frequency"
        print(out.mat)
        cat("\n")
      }
      cat("Frequencies of sites with detections, extinctions, and colonizations:\n")
      ##add matrix of frequencies
      print(x$out.freqs)
  }
}
