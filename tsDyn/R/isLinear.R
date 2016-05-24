#
#
#  Linearity Tests
#
#



#'isLinear
#'
#'Generic NLAR linearity test
#'
#'
#'@param object fitted time series model
#'@param ... arguments to and from other methods
#'@author A. F. Di Narzo
#'@keywords ts
#'@export
isLinear <- function(object, ...)
  UseMethod("isLinear")

isLinear.default <- function(object, ...)
  stop("no linearity tests available for this model")

#' @S3method isLinear lstar
isLinear.lstar <- function(object, mTh, thDelay = 0, thVar, trace=TRUE, ...)
{

# Reading function arguments  
  
  externThVar <- FALSE
  
  if(!missing(thDelay)) {
    
    if(thDelay >= object$m)
      stop(paste("thDelay too high: should be < m (=",object$m,")"))
    
    s_t <- object$xx[,thDelay+1]  # s_t is the thDelay-lagged series
    
  }
  else if(!missing(mTh)) {
    
    if(length(mTh) != object$m) 
      stop("length of 'mTh' should be equal to 'm'")
    
    s_t <- object$xx %*% mTh #threshold variable as combination of lags
    dim(s_t) <- NULL
    
  }
  else if(!missing(thVar)) {
    
    if(length(thVar) > nrow(object$xx)) {
      
      s_t <- thVar[1:nrow(object$xx)]
      
      if(trace)
        cat("Using only first", nrow(object$xx), "elements of thVar\n")
      
    }
    else 
      s_t <- thVar
    
    externThVar <- TRUE
    
  }
  else {
    
    if(trace)
      cat("Using default threshold variable: thDelay=0\n")
    
    s_t <- object$xx[,1]
  }

# Parameters read
  
  sampleSize <- length(object$yy);
  T <- NROW(object$xx);  # The number of lagged samples
  
  # Build the regressand vector
  y_t <- object$yy;
  
  # Build the regressors matrix
  if (externThVar)
    x_t <- cbind(1, object$xx)
  else
    x_t <- object$xx;
  
  # "1. Regress y_t on x_t and compute the residual sum of squares"
  regression1 <- lm(y_t ~ ., data=data.frame(x_t));
  SSR0 <- sum(regression1$residuals^2);
  
  # "2. Regress y_t (or regression1$resid) on x_t and x_t * s_t
  #      (first order) and compute the residual sum of squares"
  
  aux_data1 <- data.frame(y_t = y_t, a = x_t, b = x_t * s_t);
  aux_regression1 <- lm(y_t ~ ., data=aux_data1);
  SSR1 <- sum(aux_regression1$residuals^2);
  
  # 3. Compute the first order statistic
  n <- object$m + 1;
  m <- dim(aux_data1)[2] - n;
  F_1 <- ((SSR0 - SSR1) / m) / (SSR1 / (T - n - m));
  
  # Look up the statistic in the table, get the p-value
  lmStatTaylor1 <- pf(F_1, m, T - m - n, lower.tail = FALSE);
  
  # Regress y_t on the restrictions and compute the RSS
  aux_data3 <- data.frame(y_t = y_t, a = x_t, b = x_t * s_t,
                          c = x_t * s_t^2, d = x_t * s_t^3)
  aux_regression3 <- lm(y_t ~ ., data=aux_data3)
  SSR3 <- sum(aux_regression3$residuals^2);
  
  # Compute the third order statistic
  n <- object$m + 1;
  m <- dim(aux_data3)[2] - n;
  F_3 = ((SSR0 - SSR3) / m) / (SSR3 / (T - m - n));
  
  # Look up the statistic in the table, get the p-value
  lmStatTaylor3 <- pf(F_3, m, T - m - n, lower.tail = FALSE);
  
  # Regress y_t on the restrictions and compute the RSS
  aux_data5 <- data.frame(y_t = y_t, a = x_t, b = x_t * s_t,
                          c = x_t * s_t^2, d = x_t * s_t^3,
                          e = x_t * s_t^4, d = x_t * s_t^5)
  aux_regression5 <- lm(y_t ~ ., data=aux_data5)
  SSR5 <- sum(aux_regression5$residuals^2);
  
  # Compute the fifth order statistic
  n <- object$m + 1;
  m <- dim(aux_data5)[2] - n;
  F_5 = ((SSR0 - SSR5) / m) / (SSR5 / (T - m - n));
  
  # Look up the statistic in the table, get the p-value
  lmStatTaylor5 <- pf(F_5, m, T - m - n, lower.tail = FALSE);
  
  c(firstOrderTest = lmStatTaylor1, thirdOrderTest = lmStatTaylor3,
    fifthOrderTest = lmStatTaylor5)
  
}
