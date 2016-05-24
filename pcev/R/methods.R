# Permutation methods----

#' Permutation p-value
#' 
#' Computes a p-value using a permutation procedure.
#' 
#' @param pcevObj A pcev object of class \code{PcevClassical} or 
#'   \code{PcevBlock}
#' @param shrink Should we use a shrinkage estimate of the residual variance?
#' @param index If \code{pcevObj} is of class \code{PcevBlock}, index is a
#'   vector describing the block to which individual response variables
#'   correspond.
#' @param nperm The number of permutations to perform.
#' @param ... Extra parameters.
#' @export
permutePval <- function(pcevObj, ...) UseMethod("permutePval")

#' @describeIn permutePval
permutePval.default <- function(pcevObj, ...) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical or PcevBlock"),
       call. = FALSE)
}

#' @describeIn permutePval
permutePval.PcevClassical <- function(pcevObj, shrink, index, nperm, ...) {
  results <- estimatePcev(pcevObj, shrink)
  N <- nrow(pcevObj$Y)
  
  PCEV <- pcevObj$Y %*% results$weights
  initFit <- lm.fit(pcevObj$X, PCEV)
  df1 <- ncol(pcevObj$X) - 1
  df2 <- N - ncol(pcevObj$X)
  initFstat <- (sum((mean(PCEV) - initFit$fitted.values)^2)/df1)/(sum(initFit$residuals^2)/df2)
  initPval <- pf(initFstat, df1, df2, lower.tail = FALSE)
  
  permutationPvalues <- replicate(nperm, expr = {
    tmp <- pcevObj
    tmp$Y <- tmp$Y[sample(N), ]
    
    tmpRes <- try(estimatePcev(tmp, shrink, index), silent=TRUE)
    if(inherits(tmpRes, "try-error")) {
      return(NA)
    } else {
      tmpPCEV <- tmp$Y %*% tmpRes$weights
      
      tmpFit <- lm.fit(tmp$X, tmpPCEV)
      tmpFstat <- (sum((mean(tmpPCEV) - tmpFit$fitted.values)^2)/df1)/(sum(tmpFit$residuals^2)/df2)
      return(pf(tmpFstat, df1, df2, lower.tail = FALSE))
    }
  })
  
  pvalue <- mean(permutationPvalues < initPval)
  
  results$pvalue <- pvalue
  
  return(results)
}

#' @describeIn permutePval
permutePval.PcevBlock <- function(pcevObj, shrink, index, nperm, ...) {
  results <- estimatePcev(pcevObj, shrink, index)
  N <- nrow(pcevObj$Y)
  
  PCEV <- pcevObj$Y %*% results$weights
  initFit <- lm.fit(pcevObj$X, PCEV)
  df1 <- ncol(pcevObj$X) - 1
  df2 <- N - ncol(pcevObj$X)
  initFstat <- (sum((mean(PCEV) - initFit$fitted.values)^2)/df1)/(sum(initFit$residuals^2)/df2)
  initPval <- pf(initFstat, df1, df2, lower.tail = FALSE)
  
  permutationPvalues <- replicate(nperm, expr = {
    tmp <- pcevObj
    tmp$Y <- tmp$Y[sample(N), ]
    
    tmpRes <- try(estimatePcev(tmp, shrink, index), silent=TRUE)
    if(inherits(tmpRes, "try-error")) {
      return(NA)
    } else {
      tmpPCEV <- tmp$Y %*% tmpRes$weights
      
      tmpFit <- lm.fit(tmp$X, tmpPCEV)
      tmpFstat <- (sum((mean(tmpPCEV) - tmpFit$fitted.values)^2)/df1)/(sum(tmpFit$residuals^2)/df2)
      return(pf(tmpFstat, df1, df2, lower.tail = FALSE))
    }
  })
  
  
  pvalue <- mean(permutationPvalues < initPval, na.rm = TRUE)
  results$pvalue <- pvalue
  
  return(results)
}

###################
# Wilks methods----

#' Wilks' lambda exact test
#' 
#' Computes a p-value using the Wilk's Lambda. The null distribution of this
#' test statistic is only known in the case of a single covariate, and therefore
#' this is the only case implemented.
#' 
#' @param pcevObj A pcev object of class \code{PcevClassical} or 
#'   \code{PcevBlock}
#' @param shrink Should we use a shrinkage estimate of the residual variance?
#' @param index If \code{pcevObj} is of class \code{PcevBlock}, index is a
#'   vector describing the block to which individual response variables
#'   correspond.
#' @param ... Extra parameters.
#' @export
wilksPval <- function(pcevObj, ...) UseMethod("wilksPval")

#' @describeIn wilksPval
wilksPval.default <- function(pcevObj, ...) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical or PcevBlock"),
       call. = FALSE)
}

#' @describeIn wilksPval
wilksPval.PcevClassical <- function(pcevObj, shrink, index, ...) {
  results <- estimatePcev(pcevObj, shrink)
  N <- nrow(pcevObj$Y)
  p <- ncol(pcevObj$Y)
  d <- results$largestRoot
  wilksLambda <- (N-p-1)/p * d
  
  df1 <- p
  df2 <- N-p-1
  pvalue <- pf(wilksLambda, df1, df2, lower.tail = FALSE)
  results$pvalue <- pvalue
  
  return(results)
}

#' @describeIn wilksPval
wilksPval.PcevBlock <- function(pcevObj, shrink, index, ...) {
  stop(strwrap("Pcev is currently not implemented for
               estimation with blocks and an exact inference method"),
       call. = FALSE)
}

################################
# Roy's largest root methods----

#' Roy's largest root exact test
#' 
#' This function uses Johnstone's approximation to the null distribution of 
#' Roy's Largest Root statistic. It uses the a location-scale variant of the 
#' Tracy-Wildom distribution of order 1.
#' 
#' Note that if \code{shrink} is set to \code{TRUE}, the location-scale
#' parameters are estimated using a smaller number of permutations. distribution
#' 
#' @param pcevObj A pcev object of class \code{PcevClassical} or 
#'   \code{PcevBlock}
#' @param shrink Should we use a shrinkage estimate of the residual variance?
#' @param index If \code{pcevObj} is of class \code{PcevBlock}, index is a 
#'   vector describing the block to which individual response variables 
#'   correspond.
#' @param ... Extra parameters.
#' @export
roysPval <- function(pcevObj, ...) UseMethod("roysPval")

#' @describeIn roysPval
roysPval.default <- function(pcevObj, ...) {
  stop(strwrap("This function should be used with a Pcev object of class 
               PcevClassical"),
       call. = FALSE)
}

#' @describeIn roysPval
roysPval.PcevClassical <- function(pcevObj, shrink, index, ...) {
  
  results <- estimatePcev(pcevObj, shrink)
  n <- nrow(pcevObj$Y)
  p <- ncol(pcevObj$Y)
  q <- ncol(pcevObj$X) - 1
  d <- results$largestRoot
  # theta <- d / (1 + d)
  
  nuH <- q
  nuE <- n - q - 1
  s <- min(p, nuH)
  m <- 0.5 * (abs(p - nuH) - 1)
  N <- 0.5 * (nuE - p - 1)
  one_third <- 1/3
  
  N1 <- 2 * (s + m + N) + 1 
  gamma <- 2 * asin( sqrt( (s - 0.5)/N1 ) )
  phi <- 2*asin( sqrt( (s + 2*m + 0.5)/N1 ) )
  
  mu <- 2 * log(tan( 0.5*(phi + gamma)))
  sigma <- (16/N^2)^one_third * ( sin(phi+gamma)^2*sin(phi)*sin(gamma)) ^(-1*one_third)
  
  if(shrink) {
    # Estimate the null distribution using 
    # permutations and MLE
    null_dist <- replicate(25, expr = {
      tmp <- pcevObj
      tmp$Y <- tmp$Y[sample(N), ]
      
      tmpRes <- try(estimatePcev(tmp, shrink, index), silent=TRUE)
      if(inherits(tmpRes, "try-error")) {
        return(NA)
      } else {
        return(tmpRes$largestRoot)
      }
    })
    
    # Fit a location-scale version of TW distribution
    # Note: likelihood may throw warnings at some evaluations, 
    # which is OK
    oldw <- getOption("warn")
    options(warn = -1)
    res <- optim(c(mu, sigma), function(param) logLik(param, log(null_dist)), 
                 control = list(fnscale=-1))
    options(warn = oldw)
    
    mu1 <- res$par[1]
    sigma1 <- res$par[2]
    TW <- (log(d) - mu1)/sigma1
    
    pvalue <- RMTstat::ptw(TW, beta=1, lower.tail = FALSE, log.p = FALSE)
  } else {
    # TW <- (log(theta/(1-theta)) - mu)/sigma
    TW <- (log(d) - mu)/sigma
    
    pvalue <- RMTstat::ptw(TW, beta=1, lower.tail = FALSE, log.p = FALSE)
  }
  
  results$pvalue <- pvalue
  
  return(results)
  
}

#' @describeIn roysPval
roysPval.PcevBlock <- function(pcevObj, shrink, index, ...) {
  
  warning("There is no theoretical guarantee that this approach works.\nWe recommend comparing the resulting p-value with that obtained from a permutation procedure.")
  
  results <- estimatePcev(pcevObj, shrink, index)
  N <- nrow(pcevObj$Y)
  p <- ncol(pcevObj$Y)
  q <- ncol(pcevObj$X) - 1
  
  # Estimate largest root
  root_Vr_bd <- blockMatrixDiagonal(results$rootVr$first)
  
  fit <- lm.fit(cbind(pcevObj$X, pcevObj$Z), pcevObj$Y)
  Yfit <- fit$fitted.values
  fit_confounder <- lm.fit(cbind(rep_len(1, N), pcevObj$Z), pcevObj$Y)
  Yfit_confounder <- fit_confounder$fitted.values
  
  Vm <- crossprod(Yfit - Yfit_confounder, pcevObj$Y)
  
  mainMatrix <- root_Vr_bd %*% Vm %*% root_Vr_bd
  
  temp1 <- eigen(mainMatrix, symmetric=TRUE, only.values = TRUE)
  d <- temp1$values[1]
  
  # theta <- d / (1 + d)
  
  nuH <- q
  nuE <- N - q - 1
  s <- min(p, nuH)
  m <- 0.5 * (abs(p - nuH) - 1)
  n <- 0.5 * (nuE - p - 1)
  one_third <- 1/3
  
  N1 <- 2 * (s + m + n) + 1 
  gamma <- 2 * asin( sqrt( (s - 0.5)/N ) )
  phi <- 2*asin( sqrt( (s + 2*m + 0.5)/N ) )
  
  mu <- 2 * log(tan( 0.5*(phi + gamma)))
  sigma <- (16/N^2)^one_third * ( sin(phi+gamma)^2*sin(phi)*sin(gamma)) ^(-1*one_third)
  
  
  null_dist <- replicate(25, expr = {
    tmp <- pcevObj
    tmp$Y <- tmp$Y[sample(N), ]
    
    tmpRes <- try(estimatePcev(tmp, shrink, index), silent=TRUE)
    if(inherits(tmpRes, "try-error")) {
      return(NA)
    } else {
      root_Vr_tmp <- blockMatrixDiagonal(tmpRes$rootVr$first)
      
      fit <- lm.fit(cbind(tmp$X, tmp$Z), tmp$Y)
      Yfit <- fit$fitted.values
      fit_confounder <- lm.fit(cbind(rep_len(1, N), tmp$Z), tmp$Y)
      Yfit_confounder <- fit_confounder$fitted.values
      
      Vm_tmp <- crossprod(Yfit - Yfit_confounder, tmp$Y)
      
      mainMatrix_tmp <- root_Vr_tmp %*% Vm_tmp %*% root_Vr_tmp
      
      tmpEi <- eigen(mainMatrix_tmp, symmetric=TRUE, only.values = TRUE)
      d <- tmpEi$values[1]
      
      return(d)
    }
  })
    
    
  # Fit a location-scale version of TW distribution
  # Note: likelihood may throw warnings at some evaluations, 
  # which is OK
  oldw <- getOption("warn")
  options(warn = -1)
  res <- optim(c(mu, sigma), function(param) logLik(param, log(null_dist)), 
               control = list(fnscale=-1))
  options(warn = oldw)
  
  mu1 <- res$par[1]
  sigma1 <- res$par[2]
  TW <- (log(d) - mu1)/sigma1
  
  pvalue <- RMTstat::ptw(TW, beta=1, lower.tail = FALSE, log.p = FALSE)
  
  results$pvalue <- pvalue
  
  return(results)
}

##################
# Print method----

#' @export
print.Pcev <- function(x, ...) {
  # Provide a summary of the results
  N <- nrow(x$pcevObj$Y)
  p <- ncol(x$pcevObj$Y)
  q <- ncol(x$pcevObj$X)
  if(x$Wilks) {
    exact <- "Wilks' lambda test)"
  } else {
    exact <- "Roy's largest root test)"
  }
  
  cat("\nPrincipal component of explained variance\n")
  cat("\n", N, "observations,", p, "response variables\n")
  cat("\nEstimation method:", x$methods[1])
  cat("\nInference method:", x$methods[2])
  if(x$methods[2] == "exact") {
    cat("\n(performed using", exact)
  }
  pvalue <- x$pvalue
  if(pvalue == 0) {
    if(x$methods[2] == "permutation") {
    pvalue <- paste0("< ", 1/x$nperm)
    }
    if(x$methods[1] == "exact") {
      pvalue <- "~ 0"
    }
  }
  cat("\nP-value obtained:", pvalue, "\n")
  if(x$shrink) cat("\nShrinkage parameter rho was estimated at", x$rho, "\n")
  cat("\nVariable importance factors")
  if(p > 10) cat(" (truncated)\n") else cat("\n")
  cat(format(sort(x$VIMP, decreasing = TRUE)[1:10], digits = 3), 
      "\n\n")
  
}

########################################################
# Utility functions for null distribution estimation----
dtw_ls <- function(x, mu, sigma, beta=1, log=FALSE) {
  x1 <- (x - mu)/sigma
  # values for dtw are only available for the 
  # interval -10 to 6
  x1 <- as.numeric(x1 >= -10)*x1 - 10*as.numeric(x1 < -10)
  x1 <- as.numeric(x1 <= 6)*x1 - 6*as.numeric(x1 > 6)
  return(RMTstat::dtw(x1, beta, log)/sigma)
}

logLik <- function(param, data) {
  mu <- param[1]
  sigma <- param[2]
  data <- data[!is.na(data)]
  
  lL <- sum(log(dtw_ls(data, mu, sigma, log=FALSE)))
  
  return(lL)
}

blockMatrixDiagonal <- function(matrix_list) {  
  
  dimensions <- sapply(matrix_list, nrow)
  finalDimension <- sum(dimensions)
  finalMatrix <- matrix(0, nrow = finalDimension, ncol = finalDimension)
  index <- 1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1), index:(index+dimensions[k]-1)] <- matrix_list[[k]]
    index <- index + dimensions[k]
  }
  
  return(finalMatrix)
}
