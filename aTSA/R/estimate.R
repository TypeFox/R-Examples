#' Estimate an ARIMA Model
#' @description Estimates an ARIMA model for a univariate time series, including a sparse
#' ARIMA model.
#' @param x a univariate time series.
#' @param p the AR order, can be a positive integer or a vector with several positive 
#' integers. The default is \code{0}.
#' @param d the degree of differencing. The default is \code{0}.
#' @param q the MA order, can be a positive integer or a vector with several positive 
#' integers. The default is \code{0}.
#' @param PDQ a vector with three non-negative integers for specification of the seasonal 
#' part of the ARIMA model. The default is \code{c(0,0,0)}.
#' @param S the period of seasonal ARIMA model. The default is \code{NA}.
#' @param method fitting method. The default is \code{CSS-ML}.
#' @param intercept a logical value indicating to include the intercept in ARIMA model. The 
#' default is \code{TRUE}.
#' @param output a logical value indicating to print the results in R console. The default is
#' \code{TRUE}.
#' @param ... optional arguments to \code{\link{arima}} function.
#' 
#' @details This function is similar to the ESTIMATE statement in ARIMA procedure of SAS, 
#' except that it does not fit a transfer function model for a univariate time series. The 
#' fitting method is inherited from \code{\link{arima}} in \code{stats} package. To be 
#' specific, the pure ARIMA(p,q) is defined as 
#' \deqn{X[t] = \mu + \phi[1]*X[t-1] + ... + \phi[p]*X[p] + 
#'              e[t] - \theta[1]*e[t-1] - ... - \theta[q]*e[t-q].}
#' The \code{p} and \code{q} can be a vector for fitting a sparse ARIMA model. For example,
#' \code{p = c(1,3),q = c(1,3)} means the ARMA((1,3),(1,3)) model defined as 
#' \deqn{X[t] = \mu + \phi[1]*X[t-1] + \phi[3]*X[t-3] + e[t] 
#' - \theta[1]*e[t-1] - \theta[3]*e[t-3].} The \code{PDQ} controls the
#' order of seasonal ARIMA model, i.e., ARIMA(p,d,q)x(P,D,Q)(S), where S is the seasonal
#' period. Note that the difference operators \code{d} and D = \code{PDQ}[2] are different. 
#' The \code{d} is equivalent to \code{diff(x,differences = d)} and D is 
#' \code{diff(x,lag = D,differences = S)}, where the default seasonal period is 
#' \code{S = frequency(x)}. 
#' 
#' The residual diagnostics plots will be drawn.
#' 
#' @note Missing values are removed before the estimate. Sparse seasonal
#' ARIMA(p,d,q)x(P,D,Q)(S) model is not allowed.
#' @return A list with class "\code{estimate}" and the same results as 
#' \code{\link{arima}}. See \code{\link{arima}} for 
#' more details.
#' @seealso \code{\link{arima}},  \code{\link{identify}},  \code{\link{forecast}}
#' @author Debin Qiu
#' @references 
#' Brockwell, P. J. and Davis, R. A. (1996). Introduction to Time Series and Forecasting. 
#' Springer, New York. Sections 3.3 and 8.3.
#' 
#' @examples estimate(lh, p = 1) # AR(1) process
#' estimate(lh, p = 1, q = 1) # ARMA(1,1) process
#' estimate(lh, p = c(1,3)) # sparse AR((1,3)) process
#' 
#' # seasonal ARIMA(0,1,1)x(0,1,1)(12) model
#' estimate(USAccDeaths, p = 1, d = 1, PDQ = c(0,1,1))
#' @importFrom stats arima
#' @importFrom stats coef
#' @importFrom stats vcov
#' @importFrom stats residuals
#' @export 
estimate <- function(x, p = 0,d = 0,q = 0, PDQ = c(0,0,0), S = NA, 
                     method = c("CSS-ML", "ML", "CSS"),intercept = TRUE,
                     output = TRUE,...)
{
  DNAME <- deparse(substitute(x))
  method <- match.arg(method)
  if (NCOL(x) > 1)
    stop("'x' must be a univariate time series")
  if (any(p < 0, q < 0, p%%1 != 0, q%%1 != 0))
    stop("'p' or 'q' must be a non-negative integers")
  if (any(PDQ < 0, PDQ%%1 != 0))
    stop("'PDQ' must be a vector with three non-negative integers")
  if (any(!is.finite(x)))
    warning("missing values exist in 'x'")
  if (!is.na(S) && (S < 0 || S%%1 != 0)) 
     stop("'period' must be a positive integer")
  if (is.na(S)) S <- frequency(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 1L)
    stop("invalid length of 'x'")
  MNAME <- switch(method,"ML" = "Maximum Likelihood Estimation",
            "CSS-ML" = "Conditional-Sum-of-Squares & Maximum Likelihood Estimation",
            "CSS" = "Conditional-Sum-of-Squares Estimation")
  METHOD <- switch(method,"CSS-ML" = "CSS-ML","ML" = "ML", "CSS" = "CSS")
  P <- PDQ[1]
  D <- PDQ[2]
  Q <- PDQ[3]
  if (d > 0 || D > 0) intercept <- FALSE
  if (length(p) == 1 && length(q) == 1) {
    p.s <- 1:p
    q.s <- 1:q
    fix <- NULL
    tr.pars <- TRUE
    season <- list(order = c(P,D,Q),period = S)
  }
  else {
    p.s <- p
    q.s <- q
    p <- max(p.s)
    q <- max(q.s)
    m <- max(0,p) + max(0,q) + as.numeric(intercept)
    fix <- rep(NA,m)
    index <- c(if (p > 0) p.s else NULL,if (q > 0) q.s + p else NULL,
               if (intercept) m)
    fix[-index] <- 0
    tr.pars <- FALSE
    season <- list(order = c(0,0,0),period = NA)
    if (any(PDQ > 0))
      stop("sparse SARIMA(p,d,q)(P,D,Q)(S) model is not allowed")
  }
  l.p <- ifelse(p > 0, length(p.s), 0)
  l.q <- ifelse(q > 0, length(q.s), 0)
  npar <- l.p + l.q + P + Q + as.numeric(intercept)
  fit <- arima(x,order = c(p,d,q),seasonal = season, include.mean = intercept,
               method = METHOD,fixed = fix, transform.pars = tr.pars,...)
  coeffs <- coef(fit)
  coeffs <- coeffs[coeffs != 0]
  cov.mtrx <- vcov(fit)
  res <- residuals(fit)
  se.coeffs <- sqrt(diag(cov.mtrx))
  t.val <- coeffs/se.coeffs
  p.val <- 2*pmin(pt(t.val,n - npar), 1 - pt(t.val,n - npar))
  LAG <- c(if (p > 0) p.s, if (q > 0) q.s, if (P > 0) 1:P, if (Q > 0) 1:Q,
           if (intercept) 1)
  est.mtrx <- matrix(c(coeffs,se.coeffs,t.val,p.val,LAG),npar,5)
  if (d == 0 || D == 0) est.mtrx <- rbind(est.mtrx[npar,], est.mtrx[-npar,])
  colnames(est.mtrx) <- c("Estimate","S.E","t.value","p.value","Lag")
  rownames(est.mtrx) <- c(if(intercept) "MU",if (p > 0) paste("AR",p.s,seq = ""), 
                          if (q > 0) paste("MA",q.s,seq = ""), if (P > 0)
                            paste("SAR",1:P), if (Q > 0) paste("SMA",1:Q)) 
  sig.err <- sqrt(fit$sigma2)
  aic <- fit$aic
  sbc <- aic + (log(n) - 2)*npar  
  Is <- sqrt(1/diag(cov.mtrx))
  cor.mtrx <- Is*cov.mtrx*rep(Is,each = npar)
  dimnames(cor.mtrx) <- list(rownames(est.mtrx),rownames(est.mtrx)) 
  WNcheck <- ts.diag(fit)  
  if (output) {
    if (any(PDQ > 0)) 
      cat(paste("SARIMA(",p,",",d,",",q,")(",P,",",D,",",Q,")(",S,")", 
                " model is estimated for variable: ", DNAME,sep = ""),"\n\n")
    else
      cat(paste("ARIMA(",p,",",d,",",q,")"," model is estimated for variable: ", DNAME,
                sep = ""),"\n\n")
    cat(MNAME,"\n")
    print(est.mtrx,digits = 3)
    cat("-----","\n","n = ",n,"; ","'sigma' = ",sig.err,";",
        " AIC = ",aic,";"," SBC = ",sbc,"\n", sep = "")
    cat("------------------------------","\n")
    cat("Correlation of Parameter Estimates","\n")
    print(cor.mtrx,digits = 3)
    cat("------------------------------","\n")
    cat("Autocorrelation Check of Residuals","\n")
    print(WNcheck,digits = 3)
    cat("------------------------------","\n")
    cat("Model for variable:",DNAME,"\n")
    if (d > 0 || D > 0 ) cat(paste("Period(s) of Differencing: ",DNAME,
                                   "(",d,",",D,")",sep = ""),"\n\n")
    if (intercept) cat("Estimated mean:",coeffs[npar],"\n")
    if (p > 0) {
      arfctor <- NULL
      for (i in 1:length(p.s)) 
        arfctor <- paste(arfctor,ifelse(coeffs[i] > 0, " + ", " - "),
                         abs(round(coeffs[i],4))," B**(",p.s[i],")",sep = "")
      cat("AR factors: 1",arfctor,"\n",sep = "")
    }
    if (q > 0) {
      mafctor <- NULL
      for (j in 1:length(q.s)) 
        mafctor <- paste(mafctor,ifelse(coeffs[j + l.p] > 0, " + ", " - "),
                         abs(round(coeffs[j + l.p],4))," B**(",q.s[j],")",sep = "")
      cat("MA factors: 1",mafctor,"\n", sep = "")
    }
    if (P > 0) {
      sarfactor <- NULL
      for (k in 1:P) {
        sar.index <- k + l.p + l.q
        sarfactor <- paste(sarfactor, ifelse (coeffs[sar.index] > 0, " + ", " - "), 
                           abs(round(coeffs[sar.index],4)), " B**(",k*S,")",sep = "")
      }        
      cat("SAR factors: 1", sarfactor, seq = "","\n")
    }
    if (Q > 0) {
      smafactor <- NULL
      for (l in 1:Q) {
        sma.index <- l + l.p + l.q + P
        smafactor <- paste(smafactor, ifelse (coeffs[sma.index] > 0, " + ", " - "), 
                           abs(round(coeffs[sma.index],4)), " B**(",l*S,")",sep = "")
      }     
      cat("SMA factors: 1", smafactor, seq = "")
    }
  }
  class(fit) <- "estimate"
  estimate <- fit
}