#' Performance metrics for estimate of connectiviy matrix A.
#'
#' @description Computes various performance metrics for estimate of 
#'  connectiviy matrix A.
#'
#' @param trueA True connectivity matrix 
#' @param est Estimated connectivity matrix 
#' @param thres Value at which the point estimate should be thresholded, i.e. 
#' edges with coefficients smaller than thres are discarded. Can be a sequence 
#' of values.
#'
#' @return A data frame with the following columns: 
#' \itemize{
#'  \item \code{Threshold} Value at which point estimate \code{est} was thresholded.
#'  \item \code{SHD} Structural Hamming distance between \code{trueA} and \code{est}.
#'  \item \code{TPR.Recall} True positive rate / recall value
#'  \item \code{FPR} False positive rate
#'  \item \code{Precision} Precision value
#'  }
#'  
#' @examples 
#' # true A
#' p  <- 3
#' A <- diag(p)*0
#' A[1,2] <- 0.8
#' A[2,3] <- -0.8
#' A[3,1] <- 0.8
#' 
#' # say an estimated connectivity matrix is given by:
#' A.est <- matrix(rnorm(p*p, 1e-3, 1e-3), ncol = p)
#' diag(A.est) <- 0
#' A.est[1,2] <- 0.76
#' A.est[2,3] <- -0.68
#' A.est[3,1] <- 0.83
#'  
#' # compute metrics with threshold 0.25
#' metricsThreshold(A, A.est, thres = 0.25)
metricsThreshold <- function(trueA, est, thres = seq(0.01, 1, by = 0.01)){
  trueA <- as.matrix(trueA)
  est <- as.matrix(est)
  
  res <- matrix(0, nrow = length(thres), ncol = 6)
  for(t in 1:length(thres)){
    est[abs(est) < thres[t]] <- 0
    res[t,] <- c(thres[t], unlist(metrics(trueA, est)))
  }
  res <- as.data.frame(res)
  res <- res[,-3]
  colnames(res) <- c("Threshold", "SHD", "TPR.Recall", "FPR", "Precision")
  res
}


simple.hamming <- function(trueA, estA){
  A.adj <- trueA
  A.adj[A.adj != 0] <- 1
  
  estA.adj <- estA
  estA.adj[estA.adj != 0] <- 1
  
  sum(abs(A.adj - estA.adj))
}

my.shd <- function(m1, m2){
  m1[m1 != 0] <- 1
  m2[m2 != 0] <- 1
  p <- dim(m1)[2]
  shd <- 0
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  m1[ind] <- 0
  shd <- shd + length(ind)/2
  ind <- which(ds < 0)
  m1[ind] <- m2[ind]
  shd <- shd + length(ind)/2
  d <- abs(m1 - m2)
  shd + sum((d + t(d)) > 0)/2
}


metrics <- function(trueA, estA){
  trueA <- as.matrix(trueA)
  estA <- as.matrix(estA)
  estA[abs(estA) > 0] <- 1
  trueA[abs(trueA) > 0] <- 1
  diff <- as.vector(estA) -  as.vector(trueA)
  errs <- as.data.frame(table(diff))
  if(is.element(-1, errs$diff)) FN <- errs[errs$diff == "-1", "Freq"] else FN <- 0
  if(is.element(1, errs$diff)) FP <- errs[errs$diff == "1", "Freq"] else FP <- 0
  P <- sum(trueA)
  N <- nrow(trueA)^2 - sum(trueA)
  TP <- P - FN
  TN <- N - FP
  
  precision <- TP/(TP+FP)
  err <- (FP+FN)/(P+N)
  
  TPR <- TP/P
  FPR <- FP/N
  shd <- my.shd(trueA, estA)
  list(shd = shd, err = err, TPR.Recall = TPR, FPR = FPR, precision = precision)
}

