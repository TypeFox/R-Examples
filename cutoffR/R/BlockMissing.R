#' Simulate a missing vector with block missing pattern.
#' 
#' @param n the length of the vector
#' @param maxlen the maximum length of missing
#' @param cnst the constant used to smooth the block missing
#' @param prob the probability a single element in the vector gets missing
#' @return the same length vector with wanted block missing pattern
#' @export
#' 
#' @examples
#' 
#' # default setting
#' rev1 <- MissSimulation()
#' # with larger missing probability
#' rev2 <- MissSimulation(prob = 0.5)
#' sum(is.na(rev1))
#' sum(is.na(rev2))
#' 
#' ## Simulation block missing pattern in the Murray-Darling Basin rainfall data
#' BlockMissing <- function() {
#' complete.chunk <- data(complete.chunk)
#'  block.size <- 3 # scale for blocks when simulating the first part
#'  n.years <- c(12, 36, 48, 48)  # number of years for four simulation parts
#'  n.stations <- c(17, 17, 37, 24)  # number of stations for four simulation parts
#'  n.prob <- c(0.05, 0.005, 0.005, 0.0005) # probability vector for each simulation part
  
#'  part1.sim <- function() MissSimulation(n = 4*n.years[1], maxlen=12, cnst=12, n.prob[1])
#'  part2.sim <- function() MissSimulation(n = 12*n.years[2], maxlen=3, cnst=3, n.prob[2])
#'  part3.sim <- function() MissSimulation(n = 12*n.years[3], maxlen=3, cnst=3, n.prob[3])
#'  part4.sim <- function() MissSimulation(n = 12*n.years[4], maxlen=3, cnst=3, n.prob[4])
  
#'  p1 <- function() {    
#'    part1.mat <- matrix(0, nrow = 4*n.years[1], ncol = n.stations[1])
#'    for (j in 1:length(part1.mat[1, ])) {
#'      part1.mat[, j] <- part1.sim()
#'      part1.missing.mat <- matrix(0, nrow = 12*n.years[1], ncol = n.stations[1])
#'      # each block value should repeate three times to get the true missing matrix  
#'      part1.missing.mat[1:nrow(part1.missing.mat), ] <- part1.mat[rep(1:nrow(part1.mat), 
#'      each=block.size), ]
#'      part1.missing.mat[part1.missing.mat==1] <- NA
#'    }
#'    return(p1.miss = part1.missing.mat)
#'  }
#'  
#'  p2 <- function() {
#'    # simulate missing matrix part2
#'    part2.mat <- matrix(0, nrow=12*n.years[2], ncol=n.stations[2])
#'    for (j in 1:length(part2.mat[1, ])) {
#'      part2.mat[, j] <- part2.sim()
#'      part2.missing.mat <- part2.mat
#'      part2.missing.mat[part2.missing.mat==1] <- NA  
#'    }
#'    return(p2.miss = part2.missing.mat)
#'  }
#'  
#'  p3 <- function() {
#'    # simulate missing matrix part3
#'    part3.mat <- matrix(0, nrow=12*n.years[3], ncol=n.stations[3])
#'    for (j in 1:length(part3.mat[1, ])) {
#'      part3.mat[, j] <- part3.sim()
#'      part3.missing.mat <- part3.mat
#'      part3.missing.mat[part3.missing.mat==1] <- NA  
#'    }
#'    return(p3.miss = part3.missing.mat)
#'  }
#'  
#'  p4 <- function() {
#'    # simulate missing matrix part3
#'  part4.mat <- matrix(0, nrow=12*n.years[4], ncol=n.stations[4])
#'    for (j in 1:length(part4.mat[1, ])) {
#'      part4.mat[, j] <- part4.sim()
#'      part4.missing.mat <- part4.mat
#'      part4.missing.mat[part4.missing.mat==1] <- NA  
#'    }
#'   return(p4.missing=part4.missing.mat)
#'  }
#'  
#'  return(complete.sim = as.data.frame(cbind(rbind(p2(), p1()), cbind(p3(),p4()))) 
#'        + complete.chunk)
#'}
#' # NOTRUN
#' # bdata <- BlockMissing()
#' # HeatStruct(bdata)

MissSimulation <- function(n = 84, maxlen = 15, cnst = 15, prob = 0.03) {
  
  rec <- rep(NA, n)
  for (i in 1:length(rec)) {  
    if (i==1) {
      rec[i] <- count <- rbinom(1, 1, prob)
    }
    if (count == maxlen) {
      rec[i] <- count <- 0
    } 
    else if (count == 0) {
      rec[i] <- count <- rbinom(1, 1, prob)
    } 
    else {
      rec[i] <- rbinom(1, 1, (count + cnst)/(maxlen + cnst))
      count <- count + rec[i]
    } 
  }
    return(rec)
}
