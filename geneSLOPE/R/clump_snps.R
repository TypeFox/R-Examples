#' Clumping procedure for SLOPE
#'
#' Clumping procedure performed on SNPs, columns of matrix \code{X}, from
#' object of class \code{\link{screeningResult}},
#' which is an output of function \code{\link{screen_snps}}.
#' SNPs are clustered based on their correlations. For details see package vignette.
#'
#' @export
#' @param screenResult object of class screeningResult
#' @param rho numeric, minimal correlation between two SNPs to be assigned to one clump
#' @param verbose logical, if TRUE (default) progress bar is shown
#'
#' @return object of class \code{\link{clumpingResult}}
#'
clump_snps <- function(screenResult, rho = 0.5, verbose = TRUE){

  if(rho>=1 | rho <= 0)
    stop("Error: Rho has to be within range (0,1)")

  if(length(screenResult$y) != nrow(screenResult$X))
    stop("Error: Length of phenotype must match number of observations in matrix with snps")

  if(class(screenResult)!="screeningResult")
    stop("Error: parameter screenResult has to be of class screeningResult")

  if(verbose){
    message("Clumping procedure has started. Depending on
            size of your data this may take several minutes.")
    total = sqrt(ncol(screenResult$X))
    # create progress bar
    pb <- txtProgressBar(min = 0, max = total, style = 3)
  }

  suma = sum((screenResult$y-mean(screenResult$y))^2)
  n = length(screenResult$y) - 2
  pVals <- apply(screenResult$X, 2, function(x) pValComp(x,screenResult$y,n,suma))

  a <- order(pVals, decreasing = FALSE)
  notClumped <- rep(TRUE, length(a))
  clumps <- list()
  representatives <- list()

  i <- 1
  while(any(notClumped)){
    idx = a[i]
    if(notClumped[idx]){
      clump <- abs(apply(screenResult$X[,notClumped, drop=FALSE], 2, cor, screenResult$X[,idx]))>rho
      clumps[[i]] <- which(notClumped)[clump]
      representatives[[i]] <- idx
      notClumped[ which(notClumped)[clump] ] <- FALSE
    }
    i = i+1
    if(verbose)
      setTxtProgressBar(pb, sqrt(i))
  }
  if(verbose) close(pb)

  nullClumps <- sapply(representatives, is.null)
  representatives <- representatives[!nullClumps]
  clumps <- clumps[!nullClumps]

  result <- structure(
    list( X = screenResult$X[,unlist(representatives)],
          y = screenResult$y,
          SNPnumber = representatives,
          SNPclumps = clumps,
          X_info = screenResult$X_info,
          selectedSnpsNumbers = screenResult$selectedSnpsNumbers[unlist(representatives)],
          X_all = screenResult$X,
          numberOfSnps = screenResult$numberOfSnps,
          selectedSnpsNumbersScreening = screenResult$selectedSnpsNumbers,
          pVals = screenResult$pVals,
          pValMax = screenResult$pValMax),
    class="clumpingResult")
  return(result)
}
