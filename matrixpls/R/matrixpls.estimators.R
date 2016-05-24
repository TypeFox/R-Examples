#'@title Parameter estimation of a model matrix
#'
#'@description
#'
#'Estimates the parameters of a model matrix (\code{inner},
#'\code{reflective}, or \code{formative}).
#'
#'@details
#'
#'Parameter estimation functions estimate the parameters of a model matrix (\code{inner},
#'\code{reflective}, or \code{formative}). These functions can be used as \code{parametersInner},
#'\code{parametersReflective}, and \code{parametersFormative} arguments for 
#'\code{\link{params.separate}}.
#'
#'When two-stage least squares regression is applied with \code{estimator.tsls}, all
#'exogenous variables are used as instrumental varuables. There is currently no check of whether sufficient
#'number of instrumental variables are available.
#'
#'\code{estimator.plscLoadings} estimates the loadings by scaling the weights \code{W} with the
#'correction factor \code{c} presented by Dijkstra (2011). This produces a MINRES estimator,
#'which constraints the loadings to be proportional to the weights.
#'The PLSc code is loosely based on code contributed by Wenjing Huang and developed with the guidance
#'by Theo Dijkstra.
#'
#'\code{estimator.plscLoadings} estimates loadings with an unconstrained single factor model,
#'which requires at least three indicators per block. The loadings of 
#'single indicator factors are estimated as 1 and two indicator factors as estimated by the
#'square root of the indicator correlation.
#'
#'Providing \code{C} or \code{IC} allows for using disattenuated or otherwise
#'adjusted correlation matrices. If not provided, these matrices are calculated using \code{S} and
#'\code{W}.
#'
#'@inheritParams matrixpls-common
#'
#'@param modelMatrix A model matrix with dependent variables on rows, independent variables on colums, and
#'non-zero elements indicating relationships. Can be either \code{inner}, \code{reflective},
#'or \code{formative} matrix.
#'
#'@param ... All other arguments are either ignored or passed through to third party estimation functions.
#'
#'@return A matrix with estimated parameters.
#'
#'@name parameterEstimators
#'
NULL

#'@describeIn parameterEstimators parameter estimation with OLS regression. Can be applied to \code{inner}, \code{reflective},
#'or \code{formative} matrix.
#'@export

estimator.ols <- function(S, modelMatrix, W, ..., C = NULL, IC = NULL){
  
  covIV <- covDV <- NULL
  
  # Calculate the composite covariance matrix
  if(is.null(C)) C <- W %*% S %*% t(W)
  
  # Calculate the covariance matrix between indicators and composites
  if(is.null(IC)) IC <- W %*% S
  
  # Covariances between the independent variables
  for(m in list(S, C)){
    if(all(colnames(modelMatrix) %in% colnames(m))){
      covIV <- m[colnames(modelMatrix),colnames(modelMatrix)]
      break()
    }
  }
  
  # Covariances between IVs and DVS
  
  for(m in list(S, C, IC)){
    if(all(rownames(modelMatrix) %in% rownames(m)) &&
       all(colnames(modelMatrix) %in% colnames(m))){
      covDV <- m[rownames(modelMatrix),colnames(modelMatrix)]
      break()
    }
  }
  
  for(row in 1:nrow(modelMatrix)){
    
    independents <- which(modelMatrix[row,]!=0, useNames = FALSE)
    
    if(length(independents)==1){
      # Simple regresion is the covariance divided by the variance of the predictor
      modelMatrix[row,independents] <- covDV[row,independents]/covIV[independents,independents]
    }
    if(length(independents)>1){
      coefs <- solve(covIV[independents,independents],covDV[row,independents])
      modelMatrix[row,independents] <- coefs
    }
  }
  
  return(modelMatrix)
}

#'@describeIn parameterEstimators parameter estimation with two-stage least squares regression. For \code{inner} matrix only.
#'@export


estimator.tsls <- function(S, modelMatrix, W, ..., C){
  
  # Calculate the composite covariance matrix
  if(is.null(C)) C <- W %*% S %*% t(W)
  
  exog <- apply(modelMatrix == 0, 1, all)
  endog <- ! exog
  
  # Use all exogenous variables as instruments
  
  instruments <- which(exog) 
  for(row in 1:nrow(modelMatrix)){
    
    independents <- which(modelMatrix[row,]!=0, useNames = FALSE)
    
    
    # Check if this variable is endogenous
    
    if(length(independents) > 0 ){
      
      # Stage 1
      
      needInstruments <- which(modelMatrix[row,]!=0 & endog, useNames = FALSE)
      
      # This is a matrix that will transform the instruments into independent
      # variables. By default, all variables are instrumented by themselves
      
      stage1Model <- diag(nrow(modelMatrix))
      
      for(toBeInstrumented in needInstruments){
        # Regress the variable requiring instruments on its predictors excluding the current DV
        
        coefs1 <- solve(C[instruments, instruments],C[toBeInstrumented, instruments])
        
        stage1Model[toBeInstrumented, toBeInstrumented] <- 0
        stage1Model[toBeInstrumented, instruments] <- coefs1
      }
      
      
      # Stage 2
      
      C2 <- stage1Model %*% C %*% t(stage1Model)
      
      if(length(independents)==1){
        # Simple regresion is the covariance divided by the variance of the predictor
        modelMatrix[row,independents] <- C2[row,independents]/C2[independents,independents]
      }
      if(length(independents)>1){
        coefs2 <- solve(C2[independents,independents],C2[row,independents])
        modelMatrix[row,independents] <- coefs2
      }
      
    }
    # Continue to the next equation
  }
  
  return(modelMatrix)
}

#'@describeIn parameterEstimators parameter estimation with Dijkstra's (2011) PLSc correction for loadings.  For \code{reflective} matrix only.
#'@author Mikko Rönkkö, Wenjing Huang, Theo Dijkstra
#'
#'@references
#'
#' Huang, W. (2013). PLSe: Efficient Estimators and Tests for Partial Least Squares (Doctoral dissertation). University of California, Los Angeles.
#' 
#' Dijkstra, T. K. (2011). Consistent Partial Least Squares estimators for linear and polynomial factor models. A report of a belated, serious and not even unsuccessful attempt. Comments are invited. Retrieved from http://www.rug.nl/staff/t.k.dijkstra/consistent-pls-estimators.pdf
#'  
#'@example example/matrixpls.plsc-example.R
#'@export


estimator.plscLoadings <- function(S, modelMatrix, W,  ...){
  
  ab <- nrow(W) #number of blocks
  ai <- ncol(W) #total number of indicators
  
  L <- modelMatrix
  
  # Indicator indices based on the reflective model. Coerced to list to avoid the problem that
  # apply can return a list or a matrix depending on whether the number of indicators is equal
  # between the LVs
  
  p_refl <- apply(L, 2, function(x){list(which(x!=0))})
  
  # Calculation of the correlations between the PLS proxies, C:
  C <- W %*% S %*% t(W)	
  
  # Calculate the covariance matrix between indicators and composites
  IC <- W %*% S
  
  # Dijkstra's correction
  
  # Determination of the correction factors, based on (11) of Dijkstra, April 7, 2011.
  c2 <- rep(1,ab)
  for (i in 1:ab) {
    idx <- p_refl[[i]][[1]]
    if (length(idx) > 1) { # only for latent factors, no need to correct for the single indicator for the phantom LV
      c2[i] <- t(W[i,idx])%*%(S[idx,idx]-diag(diag(S[idx,idx])))%*%W[i,idx]
      c2[i] <- c2[i]/(t(W[i,idx])%*%(W[i,idx]%*%t(W[i,idx])-diag(diag(W[i,idx]%*%t(W[i,idx]))))%*%W[i,idx])
    }
  }
  
  # Dijkstra's formula seems to have problems with negative weights
  c2 <- abs(c2)
  
  c <- sqrt(c2)
  
  # Determination of consistent estimates of the loadings, see (13) of Dijkstra, April 7, 2011.
  
  for (i in 1:ab) {
    idx <- p_refl[[i]][[1]]
    
    if(length(idx) > 1){
      L[idx,i] <- c[i]*W[i,idx]
    }
  }
  
  attr(L,"c") <- c
  return(L)
}

#'@describeIn parameterEstimators parameter estimation with one indicator block at at time with exploratory
#'factor analysis using the \code{\link[psych]{fa}} function from the \code{psych} package. For \code{reflective} matrix only.
#'@param fm factoring method for estimating the factor loadings. Passed through to \code{\link[psych]{fa}}.
#'@export


estimator.efaLoadings <- function(S, modelMatrix, W,  ... , fm = "minres"){
  
  L <- modelMatrix
  
  # Indicator indices based on the reflective model. Coerced to list to avoid the problem that
  # apply can return a list or a matrix depending on whether the number of indicators is equal
  # between the LVs
  
  p_refl <- apply(L, 2, function(x){list(which(x!=0))})
  
  # Loop over factors and use EFA
  
  for (i in 1:ncol(L)) {
    idx <- p_refl[[i]][[1]]
    
    if(length(idx) == 1){ # Single indicator
      L[idx,i] <- 1
    }
    else if(length(idx) == 2){ # Two indicators
      L[idx,i] <- sqrt(S[idx,idx][2])
    }
    else if(length(idx) >= 3){ # Three or more indicators
      L[idx,i] <- psych::fa(S[idx,idx], fm = fm)$loadings
    }
  }
  return(L)
}

#'@describeIn parameterEstimators Estimates a maximum likelihood confirmatory factor analysis with \code{\link[lavaan]{lavaan}}.  For \code{reflective} matrix only.
#'@export

estimator.cfaLoadings <- function(S, modelMatrix, W, ...){
  
  L <- modelMatrix # Loading pattern
  
  hasReflIndicators <- which(apply(L!=0,2,any))
  # Loadings
  parTable <- data.frame(lhs = colnames(L)[col(L)[L!=0]], op = "=~",  rhs = rownames(L)[row(L)[L!=0]], stringsAsFactors = F)
  
  # Errors
  parTable <- rbind(parTable,data.frame(lhs = rownames(L)[row(L)[L!=0]], op = "~~",  rhs = rownames(L)[row(L)[L!=0]], stringsAsFactors = F))
  
  # Factor covariances
  if(length(hasReflIndicators)>1){
    a <- matrix(0,length(hasReflIndicators),length(hasReflIndicators))
    parTable <- rbind(parTable,data.frame(lhs = colnames(L)[hasReflIndicators][col(a)[lower.tri(a)]],
                                          op = "~~",  rhs = colnames(L)[hasReflIndicators][row(a)[lower.tri(a)]], stringsAsFactors = F))
  }
  
  # Factor variances
  parTable <- rbind(parTable,data.frame(lhs = colnames(L)[hasReflIndicators], op = "~~",  rhs = colnames(L)[hasReflIndicators], stringsAsFactors = F))
  
  parTable <- cbind(id = as.integer(1:nrow(parTable)), 
                    parTable,
                    user = as.integer(ifelse(1:nrow(parTable)<=sum(L!=0),1,0)),
                    group = as.integer(1),
                    free = as.integer(ifelse(1:nrow(parTable)<=(nrow(parTable)-length(hasReflIndicators)),1:nrow(parTable),0)),
                    ustart = as.integer(ifelse(1:nrow(parTable)<=(nrow(parTable)-length(hasReflIndicators)),NA,1)),
                    exo = as.integer(0),
                    label = "",
                    eq.id = as.integer(0),
                    unco = as.integer(ifelse(1:nrow(parTable)<=(nrow(parTable)-length(hasReflIndicators)),1:nrow(parTable),0)),
                    stringsAsFactors = FALSE)
  
  args <- list(model= parTable, sample.cov = S,
               sample.nobs = 100, # this does not matter, but is required by lavaan
               se="none",
               sample.cov.rescale = FALSE,
               meanstructure = FALSE)
  
  e <- list(...)
  f <- formals(lavaan::lavaan)
  include <- intersect(names(f), names(e))
  args[include] <- e[include]
  
  cfa.res <- do.call(lavaan::lavaan, args)
  
  L[L==1] <- lavaan::coef(cfa.res)[1:sum(L!=0)]
  return(L)
  
}


#'@describeIn parameterEstimators parameter estimation with PLS regression. For \code{inner} matrix only.
#'@export
#'@author
#'Mikko Rönkkö, Gaston Sanchez, Laura Trinchera, Giorgio Russolillo
#'
#'@details
#'A part of the code for \code{\link{estimator.plsreg}} is adopted from the \pkg{plspm} package, licenced
#'under GPL3.

#'@references
#'
#'Sanchez, G. (2013). \emph{PLS Path Modeling with R.} Retrieved from http://www.gastonsanchez.com/PLS Path Modeling with R.pdf
#'
estimator.plsreg <- function(S, modelMatrix, W, ..., C){
  
  # Calculate the composite covariance matrix
  if(is.null(C)) C <- W %*% S %*% t(W)
  
  for(row in 1:nrow(modelMatrix)){
    
    independents <- which(modelMatrix[row,]!=0)
    
    if(length(independents)>0){
      vars <- c(row,independents)
      coefs <- get_plsr1(S[vars,vars])
      modelMatrix[row,independents] <- coefs
    }
  }
  
  return(modelMatrix)
  
}

#
# Run PLS regression. Ported from PLSPM
#

get_plsr1 <-function(C, nc=NULL, scaled=TRUE)
{
  # ============ checking arguments ============
  p <- ncol(C)
  if (is.null(nc))
    nc <- p
  # ============ setting inputs ==============
  if (scaled) C <- stats::cov2cor(C)
  C.old <- C
  
  Ph <- matrix(NA, p, nc)# matrix of X-loadings
  Wh <- matrix(NA, p, nc)# matrix of raw-weights
  ch <- rep(NA, nc)# vector of y-loadings
  
  # ============ pls regression algorithm ==============
  
  for (h in 1:nc)
  {
    # Covariancese between the independent and dependent vars
    w.old <- C[2:nrow(C),1]
    w.new <- w.old / sqrt(sum(w.old^2)) # normalization
    
    # Covariances between the component and the variables
    cv.new <- C %*% c(0,w.new)
    
    # Covariances between the independents and the component
    p.new <- cv.new[1]
    
    # Covariance between the dependent and the component
    c.new <- cv.new[2:length(cv.new)]
    
    # Deflation
    C.old <- C.old - cv.new
    
    Ph[,h] <- p.new
    Wh[,h] <- w.new
    ch[h] <- c.new
    
    
  }
  Ws <- Wh %*% solve(t(Ph)%*%Wh)# modified weights
  Bs <- as.vector(Ws %*% ch) # std inner coeffs    
  Br <- Bs * C[1,1]/diag(C)[2:nrow(C)]   # inner coeffs
  
  return(Br)
}
