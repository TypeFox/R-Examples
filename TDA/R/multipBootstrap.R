multipBootstrap <- 
function(Y, B = 30, alpha = 0.05, parallel = FALSE, printProgress = FALSE) {

  if (!is.numeric(Y) || !is.matrix(Y)) {
    stop("Y should be a numeric matrix")
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1) {
    stop("B should be a positive integer")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1) {
    stop("alpha should be a number between 0 and 1")
  }
  if (!is.logical(parallel)) {
    stop("parallel should be logical")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }


  n <- nrow(Y)
  Nseq <- ncol(Y)
  MeanLand <- apply(Y, 2, mean)
  Gstar <- rep(NA,B)

  if (parallel) {
    boostLapply <- parallel::mclapply
  } else {
    boostLapply <- lapply
  }

  if (printProgress) {
    cat("Bootstrap: ")
  }
  
  width <- boostLapply(seq_len(B), FUN = function(i) { 
      xi <- stats::rnorm(n)
      BootLand <-
          xi * (Y - matrix(MeanLand, nrow = n, ncol = Nseq, byrow = TRUE))
      Gstar <- max(abs(apply(BootLand, 2, sum)/sqrt(n))) 
      if (printProgress) {
        cat(i," ")
      }
      return(Gstar)
    }) 
  
  if (printProgress) {
    cat("\n")
  }
  width <- unlist(width)
  width <- stats::quantile(width, 1 - alpha) / sqrt(n)
  
  UPband1 <- MeanLand + width
  LOWband1 <- MeanLand - width
  LOWband1[which(LOWband1 < 0)] <- 0  #set negative values of lower band =0
  Band <- cbind(LOWband1, UPband1)
  
  out <- list("width" = width, "mean" = MeanLand, "band" = Band)
  
  return(out)
}
