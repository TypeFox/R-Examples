#' @title Monte Carlo simulations
#'
#' @description This function generates monte carlo simultaions of sparse VAR and 
#' its estimation (at the moment only for VAR(1) processes).
#' @param N dimension of the multivariate time series.
#' @param nobs number of observations to be generated.
#' @param nMC number of Monte Carlo simulations.
#' @param rho base value for the covariance.
#' @param sparsity density of non zero entries of the VAR matrices.
#' @param penalty penalty function to use for LS estimation. Possible values are \code{"ENET"}, 
#' \code{"SCAD"} or \code{"MCP"}.
#' @param covariance type of covariance matrix to be used in the generation of the sparse VAR model.
#' @param options (TODO: complete)
#' @param method which type of distribution to use in the generation of the entries of the matrices.
#' 
#' @return a \code{nMc}x5 matrix with the results of the Monte Carlo estimation
 
#' @export
mcSimulations <- function(N, nobs = 250, nMC = 100, rho = 0.5, sparsity = 0.05, 
                          penalty = "ENET", covariance = "toeplitz", 
                          options = NULL, method = "normal") {

  results <- list()
  
  results$confusionMatrix <- matrix(0, nMC, 4)
  results$matrixNorms <- matrix(0, nMC, 6)
  pb <- utils::txtProgressBar(min = 0, max = nMC, style = 3)
    
  for (i in 1:nMC){

      s <- simulateVAR(nobs = nobs, N = N, rho = rho, sparsity = sparsity, covariance = covariance, method = method)
      rets <- s$data$series
      genA <- s$A[[1]]
      spRad <- max(Mod(eigen(genA)$values))
      
      res <- estimateVAR(data = rets, penalty = penalty, options = options)
      
      A <- res$A[[1]]
      estSpRad <- max(Mod(eigen(A)$values))
      
      L <- A
      L[L!=0] <- 1
      L[L==0] <- 0
      
      genL <- genA
      genL[genL!=0] <- 1
      genL[genL==0] <- 0
      
      results$confusionMatrix[i, 1:4] <- prop.table(table(Predicted = L, Real = genL))
      results$accuracy[i] <- 1 -sum(abs(L-genL))/N^2   # accuracy    -(1 - sum(genL)/N^2)
      results$matrixNorms[i, 1] <- abs(sum(L)/N^2 - sparsity) # sparsity
      results$matrixNorms[i, 2] <- l2norm(A-genA) / l2norm(genA)
      results$matrixNorms[i, 3] <- frobNorm(A-genA) / frobNorm(genA)
      results$matrixNorms[i, 4] <- res$mse
      results$matrixNorms[i, 5] <- spRad
      results$matrixNorms[i, 6] <- estSpRad
      utils::setTxtProgressBar(pb, i)
    }
  
  close(pb)
  
  results$confusionMatrix <- as.data.frame(results$confusionMatrix)
  colnames(results$confusionMatrix) <- c("TP", "FP", "FN", "TN")
  results$matrixNorms <- as.data.frame(results$matrixNorms)
  colnames(results$matrixNorms) <- c("sparDiff", "l2", "frob", "mse", "spRad", "estSpRad")
  
  return(results)
  
}

