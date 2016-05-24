#' Compute Covariance Maximized Lambda4
#' 
#' @description This code estimates maximized lambda4, a split-half reliability estimate.  The function splits the halves by specifying a two column list of paired inter-item covariances in descending order. It then calculates Guttman's lambda4 on every possible split-half while preserving the inter-item pairings. The function then returns a list of the Lambda4s and then takes the minimum, maximum, median, and mean of the list.  This calculation is most appropiately applied to tests with multiple factors.
#' 
#' @param x Can be either a data matrix or a covariance matrix.
#' @param method Can specify either "Hunt" or "Osburn".
#' @param missing How to handle missing values.
#' @param show.lambda4s If TRUE then the estimates for each split are included in the output.
#' @param show.splits If TRUE then a binary matrix is exported that describes the ways the items were split.
#' @param standardize When TRUE results are standardized by using the correlation matrix instead of the covariance matrix for computation.
#' 
#' 
#' @return
#' \item{estimates}{The mean, median, max, and min of the split-half reliabilities.}
#' \item{lambda4s}{A vector of maximized split-half reliabilities.}
#' \item{method}{The method chosen.  Either "Hunt" or "Osburn".}
#' \item{Analysis.Details}{Returns the number of variables and the number of split-half reliabilities.}
#' \item{Splits}{The binary indicators of the splits for the min, max, and median split-half reliability.}
#' \item{show.splits}{Logical argument selected to show the splits.}
#' \item{show.lambdas4s}{Logical argument selected to show the split-half reliabilities.}
#' 
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' 
#' @examples
#' cov.lambda4(Rosenberg, method="Hunt")
#' cov.lambda4(Rosenberg, method="Osburn")
#' @export

cov.lambda4<-function (x, method="Hunt", missing = "complete", show.lambda4s = FALSE, show.splits = FALSE, standardize=FALSE) 
{
    nvar <- dim(x)[2]
    n <- dim(x)[1]
    p <- dim(x)[2]
    
    sigma <- impute.cov(x, missing)
    
    if(standardize==TRUE){
      sigma <- cov2cor(sigma)
    }
    
    sigma.split <- data.frame(sigma)
    sigma.split <- sigma
    sigma.split[upper.tri(sigma.split, diag = TRUE)] <- -999999
    sigma0 <- diag(sigma) - sigma 
    
    xy <- matrix(ncol = 2, nrow = nvar/2)
    for (o in 1:(nvar/2)) {
        x.m <- which(sigma.split == max(sigma.split), arr.ind = TRUE)[1,]
        xy[o, 1] <- x.m[1]
        xy[o, 2] <- x.m[2]
        sigma.split[(x.m[1]), ] <- -999999
        sigma.split[ ,(x.m[1])] <- -999999
        sigma.split[ ,(x.m[2])] <- -999999
        sigma.split[(x.m[2]), ] <- -999999
    }
    
    Ahalf <- xy[, 1]
    Bhalf <- xy[, 2]
    lftout <- which(1:nvar %in% c(Ahalf, Bhalf) == FALSE)
    if (length(c(Ahalf, Bhalf)) != nvar) {
        Bhalf <- c(Bhalf, lftout)
    }
    Ani <- length(Ahalf)
    Bni <- length(Bhalf)
    
    if(method == "Hunt"){
      Acombs <- bin.combs(Ani)
      lencombs <- nrow(Acombs)
    
      t1t.temp <- (as.numeric(1:nvar %in% Ahalf) - 0.5) * 2
      t1t.splits <- t(matrix(data = rep(t1t.temp, lencombs), nrow = nvar, 
        ncol = lencombs))
    
      full <- cbind(Acombs, Acombs)
      if (Ani != Bni) {
          full <- cbind(full, rep(1, lencombs))
      }
    
      full[, c(Ahalf, Bhalf)] <- full[, 1:nvar]
    
      if (Ani != Bni) {
          covt <- which(sigma0[lftout, ] == max(sigma0[lftout, ]))
      }
    
      if (Ani != Bni) {
          full[, lftout] <- -t1t.temp[covt]
      }
    
      t1t.matrix <- (full * t1t.splits)/2 + 0.5
      t2.matrix <- (t(t1t.matrix) - 1) * -1
      onerow <- matrix(rep(1, nvar), nrow=1)
      onevector <- t(onerow)
      l4.vect <- rep(NA, lencombs)
      for (i in 1:lencombs) {
        l4.vect[i] <- (4 * (t1t.matrix[i, ] %*% sigma %*% t2.matrix[,i]))/(onerow %*% sigma) %*% onevector
      }
      
      Max <- max(l4.vect)
      Mean <- mean(l4.vect)
      Median <- median(l4.vect)
      Minimum <- min(l4.vect)
      
      sl4 <- sort(l4.vect)
      Min.Split <- t1t.matrix[which(l4.vect == sl4[1]), ]
      if(!is.null(nrow(Min.Split)))
        {Min.Split=Min.Split[1,]}
        {Min.Split=Min.Split}
      
      Median.Split <- t1t.matrix[which(l4.vect == sl4[round(lencombs/2)]), ]
      if(!is.null(nrow(Median.Split)))
        {Median.Split=Median.Split[1,]}
        {Median.Split=Median.Split}
      
      Max.Split <- t1t.matrix[which(l4.vect == sl4[lencombs]), ]
      if(!is.null(nrow(Max.Split)))
        {Max.Split=Max.Split[1,]}
        {Max.Split=Max.Split}
      
      Splits <- data.frame(Min.Split, Median.Split, Max.Split)
      
      count <- lencombs
      lambda4 <- data.frame(Mean, Max, Median, Minimum)
      
      Analysis.Details <- data.frame(nvar, count)
      
      result <- list(estimates = lambda4, 
                     lambda4s = l4.vect, 
                     method = method,
                     Analysis.Details = Analysis.Details, 
                     Splits = Splits, 
                     show.splits=show.splits, 
                     show.lambda4s=show.lambda4s)
    }
    
    if(method == "Osburn"){
      lencombs <- 1
      t1t <- matrix(rep(NA, nvar), nrow=1)
      t1t[Ahalf] <- 1
      t1t[Bhalf] <- 0
      t2 <- (t(t1t) - 1) * -1
      onerow <- matrix(rep(1, nvar), nrow=1)
      onevector <- t(onerow)      
      l4 <- (4 *(t1t %*% sigma %*% t2) )/(onerow %*% sigma) %*% onevector
      Splits=t1t
      Analysis.Details <- data.frame(nvar, 1)
      result <- list(l4 = l4, 
                     Analysis.Details = Analysis.Details, 
                     Splits = Splits, 
                     method = method,
                     show.splits = TRUE)
    }

    class(result)<-c("cov.lambda4")
    return(result)

}
