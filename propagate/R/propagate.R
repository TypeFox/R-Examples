propagate <- function(
expr, 
data, 
type = c("stat", "raw", "sim"), 
second.order = TRUE,
do.sim = TRUE, 
dist.sim = c("norm", "t"),
df.t = NULL,
use.cov = TRUE, 
nsim = 100000,
alpha = 0.05,
...
)
{            
  options(warn = -1)
  
  ## check for correct inputs
  type <- match.arg(type)  
  dist.sim = match.arg(dist.sim)
  
  ## version 1.0-4: convert function to expression
  if (is.function(expr)) {
    ARGS <- as.list(args(expr))
    ARGS <- ARGS[-length(ARGS)]
    VARS <- names(ARGS)    
    expr <- body(expr)
    class(expr) <- "expression"
    isFun <- TRUE
  } else isFun <- FALSE
  
  if (!is.expression(expr)) stop("'expr' must be an expression")
  if (nsim < 10000) stop("'nsim' should be >= 10000 !")
  
  if (!isFun) VARS <- all.vars(expr)  
  m <- match(VARS, colnames(data))
  if (any(is.na(m))) stop("Variable names of input dataframe and expression do not match!")
  if (length(unique(m)) != length(m)) stop("Some variable names are repetitive!")
  
  if (!is.logical(use.cov)) {
    if (!is.matrix(use.cov) || dim(use.cov)[1] != dim(use.cov)[2]) stop(paste("'use.cov' is not a ", max(m), "x", max(m), " matrix!", sep=""))
  }
  
  DATA <- as.matrix(data)
  EXPR <- expr
   
  ## create data matrix depending on input
  if (type == "raw" || type == "sim") {
    meanVAL <- apply(DATA, 2, function(x) mean(x, na.rm = TRUE))
    sdVAL <- apply(DATA, 2, function(x) sd(x, na.rm = TRUE))
  }
  if (type == "stat") {
    meanVAL <- DATA[1, ]
    sdVAL <- DATA[2, ] 
  }     
  if (type == "sim") {
    isSim <- TRUE 
    do.sim <- FALSE
  } else isSim <- FALSE
  
  ## in case of type = "raw", store n for expanded uncertainty by ws.expand
  ## else assign large number 
  if (type == "raw") nSAMPLE <- apply(DATA, 2, function(x) sum(!is.na(x)))
  else nSAMPLE <- rep(100000, ncol(DATA))
   
  ## create covariances depending on input
  if (type == "raw" || type == "sim") {
    ## if 'raw'/'sim' and use covariances, make complete covariance matrix
    SIGMA <- cov(DATA, use = "complete.obs")
    ## if 'raw'/'sim' and don't use covariances, keep only diagonal variances
    if (!use.cov) {
      SIGMA[upper.tri(SIGMA)] <- 0
      SIGMA[lower.tri(SIGMA)] <- 0
    }
  }    
   
  ## if 'stat', create covariance matrix with diagonal variances
  if (type == "stat") {    
    SIGMA <- diag(DATA[2, ]^2, nrow = length(VARS), ncol = length(VARS))    
    colnames(SIGMA) <- colnames(data)
    rownames(SIGMA) <- colnames(data)
  }  
  
  ## check for same names in 'data' and 'cov'
  if (is.matrix(use.cov)) {
    m <- match(colnames(use.cov), colnames(DATA))            
    if (any(is.na(m))) stop("Names of input dataframe and var-cov matrix do not match!")             
    if (length(unique(m)) != length(m)) stop("Some names of the var-cov matrix are repetitive!")             
    SIGMA <- use.cov
  }
  
  ## replace NA's in covariance matrix with 0's
  SIGMA[is.na(SIGMA)] <- 0         
    
  ## This will bring the variables in 'data' and 'expr' in the same 
  ## order as in the covariance matrix
  m1 <- match(colnames(SIGMA), colnames(DATA))
  meanVAL <- meanVAL[m1]
  m2 <- match(colnames(SIGMA), VARS)    

  ## Monte-Carlo simulation using multivariate normal or t-distributions
  ## or input is simulated data
  if (do.sim) {     
    if (dist.sim == "norm") datSIM <- mvrnorm(nsim, mu = meanVAL, Sigma = SIGMA, empirical = TRUE)
    if (dist.sim == "t") { 
      if (type == "raw") DF <- nrow(DATA) - 1 
      if (type != "raw" & is.null(df.t)) stop("Need degrees of freedom in 'df.t' when not supplying raw observations!")
      if (type != "raw" & is.numeric(df.t)) DF <- df.t 
      datSIM <- rtmvt(nsim, mean = meanVAL, sigma = SIGMA, df = DF)     
    }      
    colnames(datSIM) <- colnames(DATA)      
  } 
  if (isSim) datSIM <- DATA    
  
  ## user-supplied simulated data to be either evaluated 'vectorized'
  ## or 'row-wise'. The former is significantly faster.
  if(do.sim || isSim) {      
    
    ## use 'row-wise' method if 'vectorized' throws an error.
    resSIM <- try(eval(EXPR, as.data.frame(datSIM)), silent = TRUE) 
    if (inherits(resSIM, "try-error")) {
      print("Using 'vectorized' evaluation gave an error. Switching to 'row-wise' evaluation...")
      resSIM <- apply(datSIM, 1, function(x) eval(EXPR, envir = as.list(x)))     
    }
       
    ## alpha-based confidence interval of MC simulations
    confSIM <- quantile(resSIM, c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)  
          
    ## warning in case of single evaluated result
    if(do.sim && length(unique(resSIM)) == 1) print("Monte Carlo simulation gave unique repetitive values! Are all derivatives constants?")   
    
    allSIM <- cbind(datSIM, resSIM)
  } else resSIM <- datSIM <- confSIM <- allSIM <- NA  
  
  ## error propagation with first/second-order Taylor expansion
  ## first-order mean: eval(EXPR), first-order variance: G.S.t(G) 
  ## version 1.0-4: continue with NA's when differentiation not possible
  MEAN1 <- try(eval(EXPR, envir = as.list(meanVAL)), silent = TRUE)
  if (!is.numeric(MEAN1)) {
    print("There was an error in calculating the first-order mean")
    MEAN1 <- NA
  }  
      
  GRAD <- try(makeGrad(EXPR, m2), silent = TRUE)  
  if (!inherits(GRAD, "try-error")) evalGRAD <- try(sapply(GRAD, eval, envir = as.list(meanVAL)), silent = TRUE)
  else {
    GRAD <- try(numGrad(EXPR, as.list(meanVAL)), silent = TRUE)  
    evalGRAD <- try(as.vector(GRAD), silent = TRUE)
  }  
      
  VAR1 <- try(t(evalGRAD) %*% SIGMA %*% matrix(evalGRAD), silent = TRUE)  
  VAR1 <- try(as.numeric(VAR1), silent = TRUE)
  if (is.na(VAR1)) print("There was an error in calculating the first-order variance")
    
  ## second-order mean: firstMEAN + 0.5 * tr(H.S), 
  ## second-order variance: firstVAR + 0.5 * tr(H.S.H.S)
  if (second.order) {
    HESS <- try(makeHess(EXPR, m2), silent = TRUE)
    if (!inherits(HESS, "try-error"))  evalHESS <- try(sapply(HESS, eval, envir = as.list(meanVAL)), silent = TRUE)
    else HESS <- try(numHess(EXPR, as.list(meanVAL)), silent = TRUE)
    evalHESS <- try(matrix(evalHESS, ncol = length(meanVAL), byrow = TRUE), silent = TRUE)  
        
    valMEAN2 <- try(0.5 * tr(evalHESS %*% SIGMA), silent = TRUE)
    valVAR2 <- try(0.5 * tr(evalHESS %*% SIGMA %*% evalHESS %*% SIGMA), silent = TRUE)
        
    MEAN2 <- try(MEAN1 + valMEAN2, silent = TRUE)
    if (inherits(MEAN2, "try-error")) {
      print("There was an error in calculating the second-order mean")
      MEAN2 <- NA
    }
        
    VAR2 <- try(VAR1 + valVAR2, silent = TRUE)  
    if (inherits(VAR2, "try-error")) {
      print("There was an error in calculating the second-order variance")
      VAR2 <- NA
    }
  } else MEAN2 <- VAR2 <- HESS <- evalHESS <- NA
  
  ## total mean and variance  
  if (second.order) totalVAR <- VAR2  else totalVAR <- VAR1
  if (second.order) totalMEAN <- MEAN2  else totalMEAN <- MEAN1
  errorPROP <- sqrt(totalVAR)  
  
  ## confidence interval based on simulation since df.residual not available
  distPROP <- rnorm(nsim, totalMEAN, errorPROP)
  confPROP <- quantile(distPROP, c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)    
  
  outSIM <- c(Mean = mean(resSIM, na.rm = TRUE), sd = sd(resSIM, na.rm = TRUE), 
              Median = median(resSIM, na.rm = TRUE), MAD = mad(resSIM, na.rm = TRUE),
              confSIM[1], confSIM[2])
  
  outPROP <- c(Mean.1 = MEAN1, Mean.2 = MEAN2, sd.1 = sqrt(VAR1), sd.2 = sqrt(VAR2), 
               confPROP[1], confPROP[2])   
  
  OUT <- list(datSIM = datSIM, resSIM = resSIM, datPROP = distPROP,
              gradient = GRAD, evalGrad = evalGRAD, covMat = SIGMA,                         
              hessian = HESS, evalHess = evalHESS, prop = outPROP, sim = outSIM,
              ws = rbind(nSAMPLE, sdVAL))
  
  class(OUT) <- "propagate"
  return(OUT)                                     
}

