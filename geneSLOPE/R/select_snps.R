#' GWAS with SLOPE
#'
#' Performs GWAS with SLOPE on given snp matrix and phenotype.
#' At first clumping procedure is performed. Highly correlated
#' (that is stronger than parameter \emph{rho}) snps are clustered.
#' Then SLOPE is used on snp matrix which contains
#' one representative for each clump.
#'
#' @export
#' @param clumpingResult clumpProcedure output
#' @param fdr, numeric, False Discovery Rate for SLOPE
#' @param lambda lambda for SLOPE. See \code{\link[SLOPE]{create_lambda}}
#' @param sigma numeric, sigma for SLOPE
#' @param verbose logical, if TRUE progress bar is printed
#' @return object of class \code{\link{selectionResult}}
#'
#' @examples
#' \dontrun{
#' slope.result <- select_snps(clumping.result, fdr=0.1)
#' }
select_snps <- function(clumpingResult, fdr = 0.1, lambda="gaussian", sigma=NULL, verbose = TRUE){
  if(fdr>=1 | fdr <= 0){
    stop("FDR has to be within range (0,1)")
  }
  if(length(clumpingResult$y) != nrow(clumpingResult$X)){
    stop("Length of y must match
         number of rows in X")
  }

  lambda <- SLOPE::create_lambda(length(clumpingResult$y),
                                 clumpingResult$numberOfSnps, fdr, "gaussian")
  lambda <- lambda[1:ncol(clumpingResult$X)]
  lambda_diffs <- diff(lambda)
  if(any(lambda_diffs==0) & which.min(lambda_diffs==0)==1){
    warning("All lambdas are equal. SLOPE does not guarantee
            False Discovery Rate control")
  }

  if (is.null(sigma) && (nrow(clumpingResult$X) >= ncol(clumpingResult$X) + 30)) {
    selected = NULL
    repeat {
      selected.prev = selected
      sigma = c(sigma, estimate_noise(clumpingResult$X[, selected], clumpingResult$y))
      result = SLOPE::SLOPE(clumpingResult$X, clumpingResult$y, fdr = fdr,
                            lambda = lambda, sigma = tail(sigma, 1))
      selected = result$selected
      if (identical(selected, selected.prev))
        break
    }
    sigma = tail(sigma, 1)
  }

  slopeResult <- SLOPE::SLOPE(X = clumpingResult$X, y = clumpingResult$y,
                              fdr = fdr, lambda = lambda, sigma = sigma)

  selectedSNPs <- unlist(clumpingResult$SNPnumber)[slopeResult$selected]
  selectedSNPs <- sort(selectedSNPs)

  X_selected <- clumpingResult$X_all[,selectedSNPs]
  if(length(selectedSNPs)==0) {
    lm.fit.summary <- summary(lm(clumpingResult$y~1))
  } else {
    # refitting linear model
    lm.fit.summary <- summary(lm(clumpingResult$y~scale(X_selected)))
  }

  result <- structure(
    list( X = X_selected,
          effects = lm.fit.summary$coefficients[-1,1],
          R2 = lm.fit.summary$r.squared,
          selectedSNPs = selectedSNPs,
          selectedClumps = clumpingResult$SNPclumps[slopeResult$selected],
          lambda = lambda,
          y = clumpingResult$y,
          clumpRepresentatives = clumpingResult$SNPnumber,
          clumps = clumpingResult$SNPclumps,
          X_info = clumpingResult$X_info,
          X_clumps = clumpingResult$X,
          X_all = clumpingResult$X_all,
          selectedSnpsNumbers = clumpingResult$selectedSnpsNumbersScreening[selectedSNPs],
          clumpingRepresentativesNumbers = clumpingResult$selectedSnpsNumbers,
          screenedSNPsNumbers = clumpingResult$selectedSnpsNumbersScreening,
          numberOfSnps = clumpingResult$numberOfSnps,
          pValMax = clumpingResult$pValMax,
          fdr = fdr),
    class="selectionResult")
  return(result)
}

