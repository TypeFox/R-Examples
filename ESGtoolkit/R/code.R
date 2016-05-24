 
# Tools -------------------------------------------------------------------

# quantiles
quantileESG <- function (x, probs) 
{
  eps <- 100 * .Machine$double.eps
  if (any((p.ok <- !is.na(probs)) & (probs < -eps | probs > 
                                       1 + eps))) 
    stop("'probs' outside [0,1]")
  n <- length(x)
  if (na.p <- any(!p.ok)) {
    o.pr <- probs
    probs <- probs[p.ok]
    probs <- pmax(0, pmin(1, probs))
  }
  np <- length(probs)
  if (n > 0 && np > 0) {
    index <- 1 + (n - 1) * probs
    lo <- floor(index)
    hi <- ceiling(index)
    x <- sort(x, partial = unique(c(lo, hi)))
    qs <- x[lo]
    i <- which(index > lo)
    h <- (index - lo)[i]
    qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
  }
  else {
    qs <- rep(NA_real_, np)
  }
  qs
}
quantileESG <- cmpfun(quantileESG)

# normal distrib. with mean = 0, sd = 1 
rnormESG <- function(n, m = NULL)
{
    if (m == 1 || is.null(m))
    {
      return(as.vector(rnormESGcpp(N = n, M = 1)))
    }
    else 
    {
      return(rnormESGcpp(N = n, M = m))
    }
}
rnormESG <- cmpfun(rnormESG)

# simulation with TAG
TAG <- function(n, m) 
{
  TAGbase <- function(n)
  {
    n2 <- 10000
    sim <- rnormESG(n = n2)
    sj <- quantileESG(sim, (0:n)/n)
    sj_up <- sj[-1]
    sj_down <- sj[-(n+1)]
    out <- TAGcorecpp(sim = sim, sj_down = sj_down, 
                      sj_up = sj_up, n = n2, p = n) 
    sample(out)
  }

  if (m == 1 || is.null(m))
  {
    return(TAGbase(n))
  }
  else
  {
    return(t(replicate(m, TAGbase(n))))
  }
}
TAG <- cmpfun(TAG)

# scaling a matrix
scaleESG <- function (x, center = TRUE, scale = TRUE) 
{
#  x <- as.matrix(x)
  nc <- ncol(x)
  if (is.logical(center)) {
    if (center) {
      center <- colMeans(x, na.rm = TRUE)
      x <- sweep(x, 2L, center, check.margin = FALSE)
    }
  }
  else if (is.numeric(center) && (length(center) == nc)) 
    x <- sweep(x, 2L, center, check.margin = FALSE)
  else stop("length of 'center' must equal the number of columns of 'x'")
  if (is.logical(scale)) {
    if (scale) {
      f <- function(v) {
        v <- v[!is.na(v)]
        sqrt(sum(v^2)/max(1, length(v) - 1L))
      }
      scale <- apply(x, 2L, f)
      x <- sweep(x, 2L, scale, "/", check.margin = FALSE)
    }
  }
  else if (is.numeric(scale) && length(scale) == nc) 
    x <- sweep(x, 2L, scale, "/", check.margin = FALSE)
  else stop("length of 'scale' must equal the number of columns of 'x'")
  return(x)
}
scaleESG <- cmpfun(scaleESG)

# plot bands
bands.plot <- function(x, y.mean, ci.upper, ci.lower, col, y.goal = NULL, goal.col = "blue", ...)
{
  if (missing(x)) stop("'x' must be provided")
  if (missing(y.mean)) stop("'x' must be provided")
  if (missing(ci.upper) & missing(ci.lower)) stop("'ci.upper' and 'ci.lower' must be provided")
  
  plot(x = x, y = y.mean, type = 'l', ...)
  polygon(c(x, rev(x)), 
          c(ci.upper, rev(ci.lower)), 
          col = col, border = FALSE)
  lines(x, y.mean, lwd = 2)  
  if (!is.null(y.goal))
  {
    abline(h = y.goal, lty = 2, lwd = 2, col = goal.col)
  }
}

# add bands on a plot
bands.add <- function(x, y.mean, col, ci.upper, ci.lower)
{
  if (missing(x)) stop("'x' must be provided")
  if (missing(col)) stop("'col' must be provided")
  
  polygon(c(x, rev(x)), 
          c(ci.upper, rev(ci.lower)), 
          col = col, border = FALSE)
  lines(x, y.mean, lwd = 2)
}

# Martingale Tests and Monte Carlo convergence --------------------------------------------------

#'@title Stochastic discount factors or discounted values
#'
#'@description 
#'
#'This function provides calculation of stochastic discount factors 
#'or discounted values
#'
#'@details
#'
#'The function result is : 
#'
#'\deqn{X_t exp(-\int_0^t r_s ds)}
#'
#'where \eqn{X_t} is an asset value at a given maturity \eqn{t}, and 
#'\eqn{(r_s)_s} is the risk-free rate.
#'
#'@param r the short rate, a \code{numeric} (constant rate) or a time series object
#'
#'@param X the asset's price, a \code{numeric} (constant payoff or asset price) or a time series 
#'object
#'
#'@seealso \code{\link{esgmcprices}}, \code{\link{esgmccv}}
#'
#'@author
#'
#'Thierry Moudiki
#'
#'@export
#'
#'@examples 
#'
#'kappa <- 1.5
#'V0 <- theta <- 0.04
#'sigma_v <- 0.2
#'theta1 <- kappa*theta
#'theta2 <- kappa
#'theta3 <- sigma_v
#'
#'# OU
#'r <- simdiff(n = 10, horizon = 5, 
#'                frequency = "quart",  
#'                model = "OU", 
#'                x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = theta3)
#'
#'# Stochastic discount factors
#' esgdiscountfactor(r, 1)
esgdiscountfactor <- function(r, X)
{
  if(missing(r) || missing(X))
    stop("'r' and 'X' must be provided")
  
  length.r <- length(r)
  length.X <- length(X)
  start.r <- start(r)
  deltat.r <- deltat(r)
  start.X <- start(X)
  deltat.X <- deltat(X)
  
  if(length.r == 1 && length.X == 1) 
  {
    return(X*exp(-r)) 
  }
  
  if(length.r == 1 && length.X != 1) 
  {
    r <- ts(matrix(r, nrow(X), ncol(X)), 
            start = start.X, deltat = deltat.X)
    
    if(tsp(X)[1] > 0)
    {
      return(ts(X*exp(-apply(r, 2, cumsum)*deltat.X), 
                start = 0, 
                deltat = deltat.X))      
    }
    else
    {
      Int_r <- exp(-apply(r, 2, cumsum)*deltat.X)
      return(ts(X*rbind(rep(1, ncol(X)), 
                        Int_r[1:(nrow(X)-1), ]),                                           , 
                start = 0, 
                deltat = deltat.X))
    }
  } 
  
  if(length.r != 1 && length.X == 1)
  {
    X <- ts(matrix(X, nrow(r), ncol(r)), 
            start = start.r, deltat = deltat.r)
    
    return(ts(X*exp(-apply(r, 2, cumsum)*deltat.r),                                          
              start = 0, 
              deltat = deltat.r))
  }
  
  if(length.r != 1 && length.X != 1)
  {
    if(tsp(X)[1] > 0)
    {
      return(suppressWarnings(ts(X*window(exp(-apply(r, 2, cumsum)*deltat.r), 
                                          start = start.X, 
                                          deltat = deltat.X), 
                                 start = 0, 
                                 deltat = deltat.X)))
    }
    else
    {
      Int_r <- ts(exp(-apply(r, 2, cumsum)*deltat.r), deltat = deltat.r)
      return(suppressWarnings(ts(X*rbind(rep(1, ncol(X)), 
                                         Int_r[1:(nrow(X)-1), ]),                                           , 
                                 start = 0, 
                                 deltat = deltat.X)))
    }
  }
}
esgdiscountfactor <- cmpfun(esgdiscountfactor)

#'@title
#'
#'Estimation of discounted asset prices
#'
#'@description
#'
#'This function computes estimators (sample mean) of 
#'
#'\deqn{E[X_T exp(-\int_0^T r_s ds)]}
#'
#'where \eqn{X_T} is an asset value at given maturities \eqn{T}, and 
#'\eqn{(r_s)_s} is the risk-free rate.
#'
#'@param r a \code{numeric} or a time series object, the risk-free rate(s).
#'
#'@param X asset prices obtained with \code{\link{simdiff}}
#'
#'@param maturity the corresponding maturity (optional). If missing, all the maturities 
#'available in \code{X} are used.
#'
#'@seealso \code{\link{esgdiscountfactor}}, \code{\link{esgmccv}}
#'
#'@author
#'
#'Thierry Moudiki
#'
#'@export
#'
#'@examples
#'
#'# GBM
#'
#'r <- 0.03
#'
#'eps0 <- simshocks(n = 100, horizon = 5, frequency = "quart")
#'sim.GBM <- simdiff(n = 100, horizon = 5, frequency = "quart",   
#'                model = "GBM", 
#'                x0 = 100, theta1 = 0.03, theta2 = 0.1, 
#'                eps = eps0)
#'
#'# monte carlo prices
#'esgmcprices(r, sim.GBM)
#'
#'# monte carlo price for a given maturity
#'esgmcprices(r, sim.GBM, 2)
#'
esgmcprices <- function(r, X, maturity = NULL)
{
  if(missing(r) || missing(X))
    stop("'r' and 'X' must be provided")
  
  maturity.out <- maturity
  
  if(is.ts(X) && tsp(X)[1] > 0 && !is.null(maturity))
  {
    maturity.out <- maturity - deltat(X)
  }
  
  Y <- esgdiscountfactor(r, X)
  
  if(length(r) == 1 && length(X) == 1) 
  {
    return(Y) 
  }
  
  Z <- ts(rowMeans(Y), start = start(Y), deltat = deltat(Y))
  
  if(!is.null(maturity))
  {
    return(window(Z, start = maturity.out, end = maturity.out))
  }
  else  
  {
    return(Z)
  }
  
}
esgmcprices <- cmpfun(esgmcprices)

#'@title
#'
#'Convergence of Monte Carlo prices
#'
#'@description
#'
#'This function computes and plots confidence intervals around the estimated 
#'average price, as functions of the number of simulations. 
#'
#'@details
#'
#'Studying the convergence of the sample mean of :
#'
#'\deqn{E[X_T exp(-\int_0^T r_s ds)]}
#'
#'towards its true value.
#'
#'@param r a \code{numeric} or a time series object, the risk-free rate(s).
#'
#'@param X asset prices obtained with \code{\link{simdiff}}
#'
#'@param maturity the corresponding maturity (optional). If missing, all the maturities 
#'available in \code{X} are used.
#'
#'@param plot if \code{TRUE} (default), a plot of the convergence is displayed.
#'
#'@param ... additional parameters provided to \code{\link{matplot}}
#'
#'@return a list with estimated average prices and the confidence intervals around them.
#'
#'@author
#'
#'Thierry Moudiki
#'
#'@examples
#'
#'r <- 0.03
#'
#'set.seed(1)
#'eps0 <- simshocks(n = 100, horizon = 5, frequency = "quart")
#'sim.GBM <- simdiff(n = 100, horizon = 5, frequency = "quart",   
#'                model = "GBM", 
#'                x0 = 100, theta1 = 0.03, theta2 = 0.1, 
#'                eps = eps0)
#'
#'# monte carlo prices
#'esgmcprices(r, sim.GBM)
#'
#'# convergence to a specific price
#'(esgmccv(r, sim.GBM, 2))
#'
esgmccv <- function(r, X, maturity, plot = TRUE, ...)
{
  if(missing(r) || missing(X) || missing(maturity))
    stop("'r', 'X', and 'maturity' must be provided")
  
  bool.X <- (is.ts(X) && tsp(X)[1] > 0)
  
  if(bool.X)
    maturity <- maturity - deltat(X)
  
  Y <- esgdiscountfactor(r, X)
  Z <- window(Y, start = maturity, end = maturity)
  N <- length(Z)  
  
  avg.price <- sapply(2:N, function(x) mean(Z[1:x]))
  conf.int <- t(sapply(2:N, function(x) t.test(Z[1:x])$conf.int[1:2]))
  colnames(conf.int) <- c("lower bound", "upper bound")
  
  if (plot == TRUE)
  {
    x <- 2:N
    matplot(x, conf.int, type = 'l', xlab = "number of simulations", 
            ylab = "monte carlo estim. price", lty = c(1, 1), lwd = 2, ...)
    polygon(c(x, rev(x)), 
            c(as.vector(conf.int[, 2]), rev(as.vector(conf.int[, 1]))), 
            col = "lightyellow", border = FALSE)
    lines(x, avg.price, col = "blue")    
  }  
  invisible(list(avg.price = avg.price,
                 conf.int = conf.int))                  
}
esgmccv <- cmpfun(esgmccv)


#'@title Martingale and market consistency tests
#'
#'@description
#'
#'This function performs martingale and market consistency (t-)tests.
#'
#'@param r a \code{numeric} or a time series object, the risk-free rate(s).
#'
#'@param X a time series object, containing payoffs or projected asset values.
#' 
#'@param p0 a \code{numeric} or a vector or a univariate time series containing  
#'initial price(s) of an asset. 
#'
#'@param alpha 1 - confidence level for the test. Default value is 0.05.
#'
#'@return The function result can be just displayed. Otherwise, you can get a list 
#'by an assignation, containing (for each maturity) : 
#'\itemize{
#'\item the Student t values 
#'\item the p-values 
#'\item the estimated mean
#'of the martingale difference
#'\item  Monte Carlo prices
#'}
#'
#'@export
#'
#'@seealso \code{\link{esgplotbands}}
#'
#'@author Thierry Moudiki
#'
#'@examples
#'
#'r0 <- 0.03
#'S0 <- 100
#'
#'set.seed(10)
#'eps0 <- simshocks(n = 100, horizon = 3, frequency = "quart")
#'sim.GBM <- simdiff(n = 100, horizon = 3, frequency = "quart",   
#'                model = "GBM", 
#'                x0 = S0, theta1 = r0, theta2 = 0.1, 
#'                eps = eps0)
#' 
#'mc.test <- esgmartingaletest(r = r0, X = sim.GBM, p0 = S0, alpha = 0.05)                               
#'esgplotbands(mc.test)                
#'
esgmartingaletest <- function(r, X, p0, alpha = 0.05)
{   
  delta_X <- deltat(X)
  
  if (length(r) == 1) 
  {
    r <- ts(data = matrix(data = r, nrow = nrow(X), ncol = ncol(X)), 
            start = 0, deltat = delta_X)
  }  
  nrow.r <- nrow(r)
  ncol.r <- ncol(r)  
  delta_r <- deltat(r)  
  nb.p0 <- length(p0)
  
  if (nb.p0 == 1)
  {
    Y <- ts(matrix(data = rep(p0, nrow.r*ncol.r), 
                   nrow = nrow.r, ncol = ncol.r), 
            start = 0, deltat = delta_X)
  }
  else
  {
    Y <- ts(data = replicate(ncol.r, p0), 
            start = 0, deltat = delta_X)    
  }
  
  delta_Y <- delta_X  
  Dt <- esgdiscountfactor(r, X)
  MartingaleDiff <- Dt - Y
  
  n <- ncol(MartingaleDiff) 
  meanMartingaleDiff <- rowMeans(MartingaleDiff[-1, ])  
  sdMartingaleDiff <- apply(MartingaleDiff[-1, ], 1, sd)  
  qtStudent <- qt(p = 1 - alpha/2, df = (n - 1))  
  stat_t <- meanMartingaleDiff/(sdMartingaleDiff/sqrt(n))  
  p_value <- pt(q = abs(stat_t), df = (n - 1), lower.tail = F) + 
             pt(q = -abs(stat_t), df = (n - 1))  
  horizon <- end(r)[1]
  mc.ci <- list(t = stat_t, 
                p.value = p_value, 
                samplemean = meanMartingaleDiff, 
                conf.int = ts(cbind(c(0, meanMartingaleDiff - qtStudent * sdMartingaleDiff/sqrt(n)), 
                                 c(0, meanMartingaleDiff + qtStudent * sdMartingaleDiff/sqrt(n))),
                              start = 0, deltat = delta_Y),
                truemean = rep.int(0, dim(MartingaleDiff)[1]), 
                true_prices = Y[1:nrow(MartingaleDiff), 1], mc.prices = rowMeans(Dt))  
  start_Y <- start(Y)  
  
  mat.ci <- ts(mc.ci$conf.int, start = 0, deltat = delta_Y)  
  colnames(mat.ci) <- c("c.i lower bound", "c.i upper bound")  
  t_val_p_val <- ts(cbind(mc.ci$t, mc.ci$p.value), start = delta_Y, deltat = delta_Y)  
  colnames(t_val_p_val) <- c("t", "p-value")
  
  cat("\n")
  cat(" martingale '1=1' one Sample t-test", "\n")
  cat("\n")
  cat(" alternative hypothesis: true mean of the martingale difference is not equal to 0", 
      "\n")
  cat("\n")
  cat("df = ", n - 1)
  cat("\n")
  print(t_val_p_val)
  cat("\n")
  cat((1 - alpha) * 100, "percent confidence intervals for the mean :", "\n")
  print(mat.ci)
  invisible(mc.ci)
}


# Correlation test for shocks --------------------------------------------------------

#'@title Correlation tests for the shocks
#'
#'@description
#'
#'This function performs correlation tests for the shocks generated by \code{\link{simshocks}}.
#'
#'@param x gaussian (bivariate) shocks, with correlation, generated by \code{\link{simshocks}}.
#' 
#'@param alternative indicates the alternative hypothesis and must be one of "two.sided", 
#'"greater" or "less". 
#'
#'@param method which correlation coefficient is to be used for the test : 
#'"pearson", "kendall", or "spearman".
#'
#'@param conf.level confidence level.
#'
#'@return a list with 2 components : estimated correlation coefficients, 
#'and confidence intervals for the estimated correlations. 
#'
#'@export
#'
#'@seealso \code{\link{esgplotbands}}
#'
#'@references
#'
#'D. J. Best & D. E. Roberts (1975), Algorithm AS 89: The Upper Tail 
#'Probabilities of Spearman's rho. Applied Statistics, 24, 377-379.
#'
#'Myles Hollander & Douglas A. Wolfe (1973), Nonparametric Statistical Methods.
#' New York: John Wiley & Sons. Pages 185-194 (Kendall and Spearman tests).
#'
#'@author Thierry Moudiki + stats package 
#'
#'@examples
#'
#'nb <- 500
#'
#'s0.par1 <- simshocks(n = nb, horizon = 3, frequency = "semi",
#'family = 1, par = 0.2)
#'
#'s0.par2 <- simshocks(n = nb, horizon = 3, frequency = "semi", 
#'family = 1, par = 0.8)
#'
#'(test1 <- esgcortest(s0.par1))
#'(test2 <- esgcortest(s0.par2))
#'par(mfrow=c(2, 1))
#'esgplotbands(test1)
#'esgplotbands(test2)
#'
esgcortest <- function(x, alternative = c("two.sided", "less", "greater"),
                        method = c("pearson", "kendall", "spearman"),
                        conf.level = 0.95)
{
  y <- x[[2]]
  x <- x[[1]]
  delta_x <- deltat(x)
  delta_y <- deltat(y)
    
  if (prod(dim(x) == dim(y)) == 0) stop("We must have dim(x) == dim(y)")  
  if (delta_x != delta_y) stop("We must have deltat(x) == deltat(y)")  
  alternative <- match.arg(alternative)  
  method <- match.arg(method)  
  nbdates <- dim(x)[1]  
  return(list(cor.estimate = ts(sapply(1:nbdates, function(i) cor(x[i,], y[i,], 
                                      method = method)), 
                                start = delta_x, deltat = delta_x), 
              conf.int = ts(t(sapply(1:nbdates, function(i) cor.test(x[i,], y[i,],
                            alternative = alternative, method = method, 
                            conf.level = conf.level)$conf.int)), start = delta_x, 
                           deltat = delta_x)))
}
esgcortest <- cmpfun(esgcortest)


# Misc plots --------------------------------------------------------------


#'@title Plot time series percentiles and confidence intervals
#'
#'@description
#'
#'This function plots colored bands for time series percentiles and confidence 
#'intervals. You can use it for outputs from \code{link{simdiff}}, 
#'\code{link{esgmartingaletest}}, \code{link{esgcortest}}.
#'
#'@param x a times series object 
#'
#'@param ... additionnal (optional) parameters provided to \code{plot}
#'
#'@export
#'
#'@author Thierry Moudiki
#'
#'@seealso \code{\link{esgplotts}}
#'
#'@examples
#'
#'# Times series
#'
#'kappa <- 1.5
#'V0 <- theta <- 0.04
#'sigma <- 0.2
#'theta1 <- kappa*theta
#'theta2 <- kappa
#'theta3 <- sigma
#'x <- simdiff(n = 100, horizon = 5, 
#'frequency = "quart",  
#'model = "OU", 
#'x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = theta3)
#'
#'par(mfrow=c(2,1))
#'esgplotbands(x, xlab = "time", ylab = "values")
#'matplot(time(x), x, type = 'l', xlab = "time", ylab = "series values")
#'
#'# Martingale test
#'
#'r0 <- 0.03
#'S0 <- 100
#'sigma0 <- 0.1
#'nbScenarios <- 100
#'horizon0 <- 10
#'eps0 <- simshocks(n = nbScenarios, horizon = horizon0, frequency = "quart",
#' method = "anti")
#'sim.GBM <- simdiff(n = nbScenarios, horizon = horizon0, frequency = "quart",   
#'                  model = "GBM", 
#'                  x0 = S0, theta1 = r0, theta2 = sigma0, 
#'                  eps = eps0)
#'
#'mc.test <- esgmartingaletest(r = r0, X = sim.GBM, p0 = S0, alpha = 0.05)   
#'esgplotbands(mc.test)
#'
#'# Correlation test
#'
#'nb <- 500
#'
#'s0.par1 <- simshocks(n = nb, horizon = 3, frequency = "semi",
#'family = 1, par = 0.2)
#'
#'s0.par2 <- simshocks(n = nb, horizon = 3, frequency = "semi", 
#'family = 1, par = 0.8)
#'
#'(test1 <- esgcortest(s0.par1))
#'(test2 <- esgcortest(s0.par2))
#'par(mfrow=c(2, 1))
#'esgplotbands(test1)
#'esgplotbands(test2)
esgplotbands <- function(x, ...)
{
  if (is.ts(x))
  {    
     nrow.x <- nrow(x)      
      x0 <- x[1, 1]
     x1 <- (as.numeric(x[nrow.x, 1]))
     x2 <- (as.numeric(x[nrow.x, 2]))
     x3 <- (as.numeric(x[nrow.x, 3]))
 
     cond1 <- (x1 == x2)
     cond2 <- (x2 == x3)
     
      if(cond1 && cond2) 
      {
        x_up <- x[-c(1, nrow.x), ]
      }
     else
     {
        x_up <- x[-1, ]
     }
      
      qt.95 <- c(x0, apply(x_up, 1, function(x) quantile(x, 0.95)))
      qt.05 <- c(x0, apply(x_up, 1, function(x) quantile(x, 0.05)))
      
      x.summary <- cbind(rep(x0, 5), 
                         apply(x_up, 1, function(x) summary(x))[-3, ])
      x.ci <- cbind(rep(x0, 3), 
                    rbind(apply(x_up, 1, function(x) t.test(x)$conf.int)[1, ],
                          rowMeans(x_up),
                          apply(x_up, 1, function(x) t.test(x)$conf.int)[2, ]))
      jet.colors <- colorRampPalette( c("lightyellow", "lightgreen") )
      nbcol <- 3
      color <- jet.colors(nbcol)
      
     if(cond1 && cond2) 
     {
       abs <- as.numeric(time(x))[-(nrow.x - 1)]
     }
     else{
       abs <- as.numeric(time(x))
     }
      
      y.mean <- x.summary[3, ]
      
      bands.plot(abs, y.mean, ci.upper = x.summary[1, ], ci.lower = x.summary[5, ], 
                 col = color[1], ylim = c(min(x.summary[1,]), max(x.summary[5,])), ...)
      bands.add(abs, y.mean, col = color[2], ci.upper = qt.95, 
                ci.lower = qt.05)
      bands.add(abs, y.mean, col = color[3], ci.upper = x.summary[2, ], 
                ci.lower = x.summary[4, ])  
    }
  
  if (is.list(x))
  {
    is.cor.test <- is.numeric(try(x$cor.estimate, silent = TRUE))
    
    if (is.cor.test)
    {
      conf.int <- x$conf.int
      abs <- as.numeric(time(conf.int))
      bands.plot(abs, x$cor.estimate, ci.upper = conf.int[ , 2], ci.lower = conf.int[ , 1], 
                 col = "#C7F6B8", ylim = c(min(conf.int[ , 1]), max(conf.int[ , 2])), 
                 main = "conf. int for the correlations", xlab = "time", ylab = "conf. int.", ...)
      points(abs, x$cor.estimate, pch = 16)
    }
    
    if (!is.cor.test)
    {
      abs <- as.numeric(time(x$conf.int))
      y.mean <- rep(0, length(abs))
      
      par(mfrow = c(2, 1))
      bands.plot(abs, y.mean, ci.upper = x$conf.int[, 2], ci.lower =  x$conf.int[, 1], 
                 col = "gray80", ylim = c(min(x$conf.int[, 1]), max(x$conf.int[, 2])), 
                 xlab = "time", ylab = "conf. int.", 
                 main = "conf. int. \n for the martingale difference", ...)
      lines(abs, y.mean, col = "blue", lty = 2)
      
      plot(abs, x$true_prices, type = 'l', col = "black", 
           ylim = x$true_prices[1]*c(1 - 0.03, 1 + 0.03), 
           main = "true (black) vs \n monte carlo (blue) prices", 
           xlab = "time", ylab = "prices")
      lines(abs, x$mc.prices, col = "blue")
      points(abs, x$true_prices, col = "black", pch = 16)
      points(abs, x$mc.prices,  col = "blue",  pch = 16)
    }
  }
}



#'@title Plot time series objects
#'
#'@description This function plots outputs from \code{\link{simdiff}}.
#'
#'@details For a large number of simulations, it's preferable to use 
#'\code{\link{esgplotbands}} for a synthetic view by percentiles.
#'
#'@param x a time series object, an output from \code{\link{simdiff}}.
#'
#'@seealso \code{\link{simdiff}}, \code{\link{esgplotbands}}
#'
#'@author
#'
#'Thierry Moudiki
#'
#'@export
#'
#'@references
#'
#'H. Wickham (2009), ggplot2: elegant graphics for data analysis. Springer 
#'New York.
#'
#'@examples
#'kappa <- 1.5
#'V0 <- theta <- 0.04
#'sigma <- 0.2
#'theta1 <- kappa*theta
#'theta2 <- kappa
#'theta3 <- sigma
#'x <- simdiff(n = 10, horizon = 5, frequency = "quart",  
#'model = "OU", 
#'x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = theta3)
#'
#'esgplotts(x)
#'
esgplotts <- function(x)
{ 
  x0 <- start(x)[1]
  dt_x <- deltat(x)
  tt <- cumsum(c(x0, rep.int(dt_x, dim(x)[1]-1)))
  meltdf <- melt(cbind(tt, as.data.frame(x)), id="tt")
  group <- rep(LETTERS[1:dim(x)[2]], each = dim(x)[1])
  qplot(tt, meltdf$value, data = meltdf, geom = "line", colour=group) + 
    xlab("Maturity") + ylab("Values") + 
    theme(legend.position = "none")
}
esgplotts <- cmpfun(esgplotts)



#'@title Visualize the dependence between 2 gaussian shocks
#'
#'@description
#'
#'This function helps you in visualizing the dependence between 2 gaussian
#' shocks. 
#'
#'@param x an output from \code{\link{simshocks}}, a list with 2 components.
#'
#'@param y an output from \code{\link{simshocks}}, a list with 2 components 
#'(Optional). 
#'
#'@export
#'
#'@seealso \code{\link{simshocks}}
#'
#'@author Thierry Moudiki + some nice blogs :) 
#'
#'@references
#'
#'H. Wickham (2009), ggplot2: elegant graphics for data analysis. Springer 
#'New York.
#'
#'@examples
#'
#'# Number of risk factors
#'d <- 2
#'
#'# Number of possible combinations of the risk factors
#'dd <- d*(d-1)/2
#'
#'# Family : Gaussian copula 
#'fam1 <- rep(1,dd)
#'# Correlation coefficients between the risk factors (d*(d-1)/2)
#'par0.1 <- 0.1
#'par0.2 <- -0.9
#'
#'# Family : Rotated Clayton (180 degrees)
#'fam2 <- 13
#'par0.3 <- 2
#'
#'# Family : Rotated Clayton (90 degrees)
#'fam3 <- 23
#'par0.4 <- -2
#'
#'# number of simulations
#'nb <- 500
#'
#'# Simulation of shocks for the d risk factors
#'s0.par1 <- simshocks(n = nb, horizon = 4, 
#'family = fam1, par = par0.1)
#'
#'s0.par2 <- simshocks(n = nb, horizon = 4, 
#'family = fam1, par = par0.2)
#'
#'s0.par3 <- simshocks(n = nb, horizon = 4, 
#'family = fam2, par = par0.3)
#'
#'s0.par4 <- simshocks(n = nb, horizon = 4, 
#'family = fam3, par = par0.4)
#'
#'\dontrun{
#'esgplotshocks(s0.par1, s0.par2)
#'esgplotshocks(s0.par2, s0.par3)
#'esgplotshocks(s0.par2, s0.par4)
#'esgplotshocks(s0.par1, s0.par4)} 
#'
esgplotshocks <-  function(x, y = NULL)
{
  x <- matrix(unlist(x), ncol = 2)
  if (!is.null(y))
  {
    y <- matrix(unlist(y), ncol = 2)
    if (length(x) != length(y)) stop("'x' and 'y' must have the same dimensions")
    xvar <- c(x[,1], y[,1])
    yvar <- c(x[,2], y[,2])
    nb  <- dim(x)[1]
    zvar <- as.factor(c(rep("x", nb), rep("y", nb)))
    xy <- data.frame(xvar, yvar, zvar)
  }
  else 
  {
    xvar <- x[,1]
    yvar <- x[,2]
    nb  <- dim(x)[1]
    zvar <- as.factor(rep("x", nb))
    xy <- data.frame(xvar, yvar, zvar)
  }
  
  #placeholder plot - prints nothing at all
  empty <- ggplot()+geom_point(aes(1,1), colour="white") +
    theme(                              
      plot.background = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
  
  #scatterplot of x and y variables
  scatter <- ggplot(xy,aes(xvar, yvar)) + 
    geom_point(aes(color=zvar)) + 
    scale_color_manual(values = c("blue", "red")) + 
    theme(legend.position=c(1,1),legend.justification=c(1,1)) 
  
  #marginal density of x - plot on top
  plot_top <- ggplot(xy, aes(xvar, fill=zvar)) + 
    geom_density(alpha=.5) + 
    scale_fill_manual(values = c("blue", "red")) + 
    theme(legend.position = "none")
  
  #marginal density of y - plot on the right
  plot_right <- ggplot(xy, aes(yvar, fill=zvar)) + 
    geom_density(alpha=.5) + 
    coord_flip() + 
    scale_fill_manual(values = c("blue", "red")) + 
    theme(legend.position = "none") 
  
  #arrange the plots together, with appropriate height and width for each row and column
  grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
}
esgplotshocks  <- cmpfun(esgplotshocks)



# simulation of gaussian shocks --------------------------------------------

#'@title 
#'
#'Underlying gaussian shocks for risk factors' simulation.
#'
#'@description 
#'
#'This function makes simulations of correlated or dependent gaussian shocks for risk factors.
#'
#'@details The function shall be used along with \code{\link{simdiff}}, in order to embed  
#'correlated or dependent random gaussian shocks into simulated diffusions. 
#'\code{\link{esgplotshocks}} can help in visualizing the type of dependence 
#'between the shocks. 
#'
#'@param n number of independent observations for each risk factor.
#'
#'@param horizon horizon of projection.
#'
#'@param frequency either "annual", "semi-annual", "quarterly", "monthly", 
#'"weekly", or "daily" (1, 1/2, 1/4, 1/12, 1/52, 1/252).
#'
#'@param method either classic monte carlo, antithetic variates, moment matching, 
#'hybrid antithetic variates + moment matching or "TAG" (see the 4th reference for 
#'the latter). 
#'
#'@param family the same as \code{\link{CDVineSim}} from package \code{CDVine}. 
#'A d*(d-1)/2 integer vector of C-/D-vine pair-copula families with values 
#' 0 = independence copula, 
#' 1 = Gaussian copula, 
#' 2 = Student t copula (t-copula), 
#' 3 = Clayton copula, 
#' 4 = Gumbel copula, 
#' 5 = Frank copula, 
#' 6 = Joe copula, 
#' 7 = BB1 copula, 
#' 8 = BB6 copula, 
#' 9 = BB7 copula, 
#' 10 = BB8 copula, 
#' 13 = rotated Clayton copula (180 degrees; "survival Clayton"), 
#' 14 = rotated Gumbel copula (180 degrees; "survival Gumbel"), 
#' 16 = rotated Joe copula (180 degrees; "survival Joe"), 
#' 17 = rotated BB1 copula (180 degrees; "survival BB1"),
#' 18 = rotated BB6 copula (180 degrees; "survival BB6"),
#' 19 = rotated BB7 copula (180 degrees; "survival BB7"),
#' 20 = rotated BB8 copula (180 degrees; "survival BB8"),
#' 23 = rotated Clayton copula (90 degrees), 
#' 24 = rotated Gumbel copula (90 degrees),
#' 26 = rotated Joe copula (90 degrees), 
#' 27 = rotated BB1 copula (90 degrees), 
#' 28 = rotated BB6 copula (90 degrees), 
#' 29 = rotated BB7 copula (90 degrees), 
#' 30 = rotated BB8 copula (90 degrees), 
#' 33 = rotated Clayton copula (270 degrees), 
#' 34 = rotated Gumbel copula (270 degrees), 
#' 36 = rotated Joe copula (270 degrees), 
#' 37 = rotated BB1 copula (270 degrees), 
#' 38 = rotated BB6 copula (270 degrees), 
#' 39 = rotated BB7 copula (270 degrees), 
#' 40 = rotated BB8 copula (270 degrees)  
#'
#'@param par the same as \code{\link{CDVineSim}} from package \code{CDVine}. 
#'A d*(d-1)/2 vector of pair-copula parameters.
#'
#'@param par2 the same as \code{\link{CDVineSim}} from package \code{CDVine}. 
#'A d*(d-1)/2 vector of second parameters for pair-copula families with two 
#'parameters (t, BB1, BB6, BB7, BB8; no default).
#'
#'@param type type of the vine model:
#' 1 : C-vine
#' 2 : D-vine
#'
#'@return
#'If \code{family} and \code{par} are not provided, a univariate time 
#'series object with simulated gaussian shocks for one risk factor. Otherwise, 
#'a list of time series objects, containing gaussian shocks for each risk factor. 
#'
#'@author Thierry Moudiki
#'
#'@references
#'
#'Brechmann, E., Schepsmeier, U. (2013). Modeling Dependence with C-
#'and D-Vine Copulas: The R Package CDVine. Journal of Statistical Software,
#'52(3), 1-27. URL \url{http://www.jstatsoft.org/v52/i03/}.
#'
#'Genz, A. Bretz, F., Miwa, T. Mi, X., Leisch, F., Scheipl, F., Hothorn, T. (2013).
#' mvtnorm: Multivariate Normal and t Distributions. R package version 0.9-9996.
#'
#'Genz, A. Bretz, F. (2009), Computation of Multivariate Normal and t Probabilities. 
#'Lecture Notes in Statistics, Vol. 195., Springer-Verlag, Heidelberg. ISBN 978-3-642-01688-2.
#'
#'Nteukam T, O., & Planchet, F. (2012). Stochastic evaluation of life insurance 
#'contracts: Model point on asset trajectories and measurement of the error 
#'related to aggregation. Insurance: Mathematics and Economics, 51(3), 624-631.
#'URL \url{http://www.ressources-actuarielles.net/EXT/ISFA/1226.nsf/0/ab539dcebcc4e77ac12576c6004afa67/$FILE/Article_US_v1.5.pdf}
#'
#'@export
#'
#'@seealso \code{\link{simdiff}}, \code{\link{esgplotshocks}}
#'
#'@examples
#'
#'# Number of risk factors
#'d <- 6
#'
#'# Number of possible combinations of the risk factors
#'dd <- d*(d-1)/2
#'
#'# Family : Gaussian copula for all
#'fam1 <- rep(1,dd)
#'
#'# Correlation coefficients between the risk factors (d*(d-1)/2)
#'par1 <- c(0.2,0.69,0.73,0.22,-0.09,0.51,0.32,0.01,0.82,0.01,
#'         -0.2,-0.32,-0.19,-0.17,-0.06)
#'
#'                  
#'# Simulation of shocks for the 6 risk factors
#'simshocks(n = 10, horizon = 5, family = fam1, par = par1)
#'
#'
#'# Simulation of shocks for the 6 risk factors
#'# on a quarterly basis
#'simshocks(n = 10, frequency = "quarterly", horizon = 2, family = fam1, 
#'par = par1)
#'
#'
#'# Simulation of shocks for the 6 risk factors simulation
#'# on a quarterly basis, with antithetic variates and moment matching. 
#'s0 <- simshocks(n = 10, method = "hyb", horizon = 4, 
#'family = fam1, par = par1)
#'
#'  
#' s0[[2]]
#' colMeans(s0[[1]])
#' colMeans(s0[[5]])
#' apply(s0[[3]], 2, sd)
#' apply(s0[[4]], 2, sd)
#'
simshocks <- function(n, horizon, 
                      frequency = c("annual", "semi-annual", 
                                    "quarterly", "monthly", 
                                    "weekly", "daily"), 
                      method = c("classic", "antithetic", 
                                 "mm", "hybridantimm", "TAG"), 
                      family = NULL, par = NULL, par2 = NULL, type = c("CVine", "DVine"))
{
  if (floor(n) != n) stop("'n' must be an integer")
  if (n <= 1) stop("'n' must be > 1")
  if (floor(horizon) != horizon) stop("'horizon' must be an integer")
  frequency <- match.arg(frequency)
  delta <- switch(frequency, 
                  "annual" = 1, 
                  "semi-annual" = 0.5, 
                  "quarterly" = 0.25, 
                  "monthly" = 1/12,
                  "weekly" = 1/52,
                  "daily" = 1/252)
  method <- match.arg(method)
  m <- horizon/delta  
  type <- match.arg(type)
  
  if (is.null(family) && is.null(par))
  {
    if (method == "classic")
    {
      return(ts(data = rnormESGcpp(N = n, M = m), 
                start = delta, end = horizon, deltat = delta))  
    }
    
    if (method == "antithetic")
    {
      half.n <- n/2
      if(floor(half.n) != half.n) stop("'n' must be divisible by 2")        
      temp <- rnormESGcpp(N = half.n, M = m)        
      return(ts(data = cbind(temp, -temp), 
                start = delta, end = horizon, deltat = delta))
    }
    
    if (method == "mm") 
    {
      return(ts(data = scaleESG(rnormESGcpp(N = n, M = m)), 
                start = delta, end = horizon, deltat = delta))  
    }
    
    if (method == "hybridantimm")
    {
      half.n <- n/2
      if(floor(half.n) != half.n) stop("'n' must be divisible by 2")
      temp <- rnormESGcpp(N = half.n, M = m)
      return(ts(data = scaleESG(cbind(temp, -temp)), 
                start = delta, end = horizon, deltat = delta))
    }
    
    if (method == "TAG")
    {
      return(ts(data = TAG(n, m), 
                start = delta, end = horizon, deltat = delta))        
    }    
  } 
  else # !is.null(family) && !is.null(par)
  { 
    nb.sim <- n*m
    
    if (method == "classic")
    {
      shocks.sim <- qnorm(CDVineSim(N = nb.sim, 
                                    family, par, par2, type)) 
      d <- dim(shocks.sim)[2]  
      return(lapply(1:d, function(i) 
        ts(data = matrix(shocks.sim[ , i], nrow = m, ncol = n), 
           start = delta, end = horizon, deltat = delta)))        
    }
    
    if (method == "antithetic")
    {
      half.nb.sim <- nb.sim/2
      if(floor(half.nb.sim) != half.nb.sim) stop("'n' must be divisible by 2")      
      temp <- qnorm(CDVineSim(N = half.nb.sim, 
                              family, par, par2, type))         
      shocks.sim <- rbind(temp, -temp)
      d <- dim(shocks.sim)[2]  
      return(lapply(1:d, function(i) 
        ts(data = matrix(shocks.sim[ , i], nrow = m, ncol = n), 
           start = delta, end = horizon, deltat = delta)))
    }
    
    if (method == "mm") 
    {
      shocks.sim <- qnorm(CDVineSim(N = nb.sim, 
                                    family, par, par2, type)) 
      d <- dim(shocks.sim)[2]  
      return(lapply(1:d, function(i) 
        ts(data = scaleESG(matrix(shocks.sim[ , i], nrow = m, ncol = n)), 
           start = delta, end = horizon, deltat = delta))) 
    }
    
    if (method == "hybridantimm")
    {
      half.nb.sim <- nb.sim/2
      if(floor(half.nb.sim) != half.nb.sim) stop("'n' must be divisible by 2")
      temp <- qnorm(CDVineSim(N = half.nb.sim, 
                              family, par, par2, type)) 
      shocks.sim <- rbind(temp, -temp)
      d <- dim(shocks.sim)[2]  
      return(lapply(1:d, function(i) 
        ts(data = scaleESG(matrix(shocks.sim[ , i], nrow = m, ncol = n)), 
           start = delta, end = horizon, deltat = delta)))  
    }
    
    if (method == "TAG")
    {
      if (!is.null(family) && family > 0) warning("for method == 'TAG', only the independence copula is implemented")   
      return(lapply(1:d, function(i) 
        ts(data = TAG(n = n, m = m), 
           start = delta, end = horizon, deltat = delta)))            
    }    
  }
}
simshocks <- cmpfun(simshocks)


# simulation of risk factors ----------------------------------------------

#'@title 
#'
#'Simulation of diffusion processes. 
#'
#'@description 
#'
#'This function makes simulations of diffusion processes, that are building 
#'blocks for various risk factors' models. 
#'
#'@param n number of independent observations.
#'
#'@param horizon horizon of projection.
#'
#'@param frequency either "annual", "semi-annual", "quarterly", "monthly", 
#'"weekly", or "daily" (1, 1/2, 1/4, 1/12, 1/52, 1/252).
#'
#'@param model either Geometric Brownian motion-like (\code{"GBM"}), 
#'Cox-Ingersoll-Ross (\code{"CIR"}), or Ornstein-Uhlenbeck 
#'(\code{"OU"}).
#'
#' GBM-like (GBM, Merton, Kou, Heston, Bates)
#'
#'\deqn{dX_t = \theta_1(t) X_t dt + \theta_2(t) X_t dW_t + X_t JdN_t}
#'  
#'  CIR
#'
#'\deqn{dX_t = (\theta_1 - \theta_2 X_t) dt + \theta_3\sqrt(X_t) dW_t}
#'  
#'  Ornstein-Uhlenbeck
#'\deqn{dX_t = (\theta_1 - \theta_2 X_t)dt + \theta_3 dW_t}
#'  
#'Where \eqn{(W_t)_t} is a standard brownian motion :
#'  
#'  \deqn{dW_t ~~ \epsilon \sqrt(dt)}
#'  
#'and \deqn{\epsilon ~~ N(0, 1)}
#'
#'The \eqn{\epsilon} is a gaussian increment that 
#'can be an output from \code{\link{simshocks}}. 
#'
#'For 'GBM-like', \eqn{\theta_1} and \eqn{\theta_2} can be held constant, and the
#' jumps part \eqn{JdN_t} is optional. In case the jumps are used, they arise 
#' following a Poisson process \eqn{(N_t)}, with intensities \eqn{J} drawn either 
#'from lognormal or asymmetric double-exponential distribution. 
#'
#'@param x0 starting value of the process. 
#'
#'@param theta1 a \code{numeric} for \code{model = "GBM"}, \code{model = "CIR"},
#'\code{model = "OU"}. Can also be a time series object (an output from 
#'\code{simdiff} with the same number of scenarios, horizon and frequency) for 
#'\code{model = "GBM"}, and time-varying parameters. 
#'
#'@param theta2 a \code{numeric} for \code{model = "GBM"}, \code{model = "CIR"},
#'\code{model = "OU"}. Can also be a time series object (an output from 
#'\code{simdiff} with the same number of scenarios, horizon and frequency) 
#'for \code{model = "GBM"}, and time-varying parameters.
#'
#'@param theta3 a \code{numeric}, volatility for \code{model = "CIR"} and 
#'\code{model = "OU"}.
#'
#'@param lambda intensity of the Poisson process counting the jumps. Optional.
#'
#'@param mu.z mean parameter for the lognormal jumps size. Optional.
#'
#'@param sigma.z standard deviation parameter for the lognormal jumps size. 
#'Optional.
#'
#'@param p probability of positive jumps. Must belong to ]0, 1[. Optional.
#'
#'@param eta_up mean of positive jumps in Kou's model. Must belong to 
#']0, 1[. Optional.
#'
#'@param eta_down mean of negative jumps. Must belong to ]0, 1[. Optional.
#'
#'@param eps gaussian shocks. If not provided, independent shocks are 
#'generated internally by the function. Otherwise, for custom shocks, 
#'must be an output from \code{\link{simshocks}}. 
#'
#'@return a time series object. 
#'
#'@seealso \code{\link{simshocks}}, \code{\link{esgplotts}}
#'
#'@author Thierry Moudiki
#'
#'@references
#'
#'Black, F., Scholes, M.S. (1973) The pricing of options and corporate liabilities, 
#'Journal of Political Economy, 81, 637-654.
#'
#'Cox, J.C., Ingersoll, J.E., Ross, S.A. (1985) A theory of the term structure of
#' interest rates, Econometrica, 53, 385-408.
#'
#'Iacus, S. M. (2009). Simulation and inference for stochastic differential 
#'equations: with R examples (Vol. 1). Springer.
#'
#'Glasserman, P. (2004). Monte Carlo methods in financial engineering (Vol. 53). 
#'Springer.
#'
#'Kou S, (2002), A jump diffusion model for option pricing, Management Sci-
#'ence Vol. 48, 1086-1101.
#'
#'Merton, R. C. (1976). Option pricing when underlying stock returns are 
#'discontinuous. Journal of financial economics, 3(1), 125-144.
#'
#'Uhlenbeck, G. E., Ornstein, L. S. (1930) On the theory of Brownian motion, 
#'Phys. Rev., 36, 823-841.
#'
#'Vasicek, O. (1977) An Equilibrium Characterization of the Term Structure, 
#'Journal of Financial Economics, 5, 177-188.
#'
#'@export
#'
#'@examples
#'
#'kappa <- 1.5
#'V0 <- theta <- 0.04
#'sigma_v <- 0.2
#'theta1 <- kappa*theta
#'theta2 <- kappa
#'theta3 <- sigma_v
#'
#'# OU
#'
#'sim.OU <- simdiff(n = 10, horizon = 5, 
#'                frequency = "quart",  
#'                model = "OU", 
#'                x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = theta3)
#'head(sim.OU)
#'par(mfrow=c(2,1))
#'esgplotbands(sim.OU, xlab = "time", ylab = "values", main = "with esgplotbands")                
#'matplot(time(sim.OU), sim.OU, type = 'l', main = "with matplot")                
#'
#'
#'# OU with simulated shocks (check the dimensions)
#'
#'eps0 <- simshocks(n = 50, horizon = 5, frequency = "quart", method = "anti")
#'sim.OU <- simdiff(n = 50, horizon = 5, frequency = "quart",   
#'                model = "OU", 
#'                x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = theta3, 
#'                eps = eps0)
#'par(mfrow=c(2,1))
#'esgplotbands(sim.OU, xlab = "time", ylab = "values", main = "with esgplotbands")                
#'matplot(time(sim.OU), sim.OU, type = 'l', main = "with matplot")
#'# a different plot
#'esgplotts(sim.OU)
#'
#'
#'# CIR
#'
#'sim.CIR <- simdiff(n = 50, horizon = 5, 
#'                frequency = "quart",  
#'                model = "CIR", 
#'                x0 = V0, theta1 = theta1, theta2 = theta2, theta3 = 0.05)
#'esgplotbands(sim.CIR, xlab = "time", ylab = "values", main = "with esgplotbands")                  
#'matplot(time(sim.CIR), sim.CIR, type = 'l', main = "with matplot")
#'
#'
#'
#'# GBM
#'
#'eps0 <- simshocks(n = 100, horizon = 5, frequency = "quart")
#'sim.GBM <- simdiff(n = 100, horizon = 5, frequency = "quart",   
#'                model = "GBM", 
#'                x0 = 100, theta1 = 0.03, theta2 = 0.1, 
#'                eps = eps0)
#'esgplotbands(sim.GBM, xlab = "time", ylab = "values", main = "with esgplotbands")                
#'matplot(time(sim.GBM), sim.GBM, type = 'l', main = "with matplot")
#'
#'
#'eps0 <- simshocks(n = 100, horizon = 5, frequency = "quart")
#'sim.GBM <- simdiff(n = 100, horizon = 5, frequency = "quart",   
#'                model = "GBM", 
#'                x0 = 100, theta1 = 0.03, theta2 = 0.1, 
#'                eps = eps0)                
#'esgplotbands(sim.GBM, xlab = "time", ylab = "values", main = "with esgplotbands")                
#'matplot(time(sim.GBM), sim.GBM, type = 'l', main = "with matplot")
simdiff <- function(n, horizon, 
                    frequency = c("annual", "semi-annual", 
                                  "quarterly", "monthly", 
                                  "weekly", "daily"), 
                    model = c("GBM", "CIR", "OU"), 
                    x0, theta1 = NULL, theta2 = NULL, theta3 = NULL,
                    lambda = NULL, mu.z = NULL, sigma.z = NULL, 
                    p = NULL, eta_up = NULL, eta_down = NULL,
                    eps = NULL)
{
  if (floor(n) != n) stop("'n' must be an integer")
  if (n <= 1) stop("'n' must be > 1")
  if (floor(horizon) != horizon) stop("'horizon' must be an integer")  
  frequency <- match.arg(frequency)
  delta <- switch(frequency, 
                  "annual" = 1, 
                  "semi-annual" = 0.5, 
                  "quarterly" = 0.25, 
                  "monthly" = 1/12,
                  "weekly" = 1/52,
                  "daily" = 1/252)
  m <- horizon/delta 
  nbdates <- m + 1
  model <- match.arg(model)
  
  if (is.null(eps))
  { 
    eps <- rnormESGcpp(N = n, M = m)
  } 
  else
  {
    if(!all.equal(dim(eps), c(m, n))) 
      stop(paste("dimensions required for 'eps' : nrow (observation dates) = "
                 , m, "\n and ncol (number of scenarios) = ", n, ".\n Use 'simshocks' function
                 with same number of observations, horizon and frequency to obtain 'eps'"))
  }
  
  if(model == "GBM")
  {
    if (!is.null(theta3)) 
      warning("unused parameter : 'theta3'")
    
    if(length(theta1) == 1) 
    {
      theta1 <- matrix(data = theta1, nrow = nbdates, ncol = n)
    }
    else
    {
      if(prod(dim(theta1) == c(nbdates, n)) == 0) 
        stop(paste("dimensions required for 'theta1' : a numeric or nrow (observation dates) = "
                   , nbdates, "\n and ncol (number of scenarios) = ", n))
    }
    
    if(length(theta2) == 1) 
    {
      theta2 <- matrix(data = theta2, nrow = nbdates, ncol = n)
    }
    else
    {
      if(prod(dim(theta2) == c(nbdates, n)) == 0) 
        stop(paste("dimensions required for 'theta2' : a numeric or nrow (observation dates) = "
                   , nbdates, "\n and ncol (number of scenarios) = ", n))
    }   
  }
  
  if (model == "CIR")
  {
    if (is.null(theta1) || is.null(theta2) || is.null(theta3))
      stop("'theta1', 'theta2' and 'theta3' must be provided")
    
    if(length(theta1) != 1 || length(theta2) != 1 || length(theta3) != 1) 
      stop(paste("length required for 'theta1', 'theta2' and 'theta3' : 1"))
    
    # rCIRESGcpp(const int N, const int horizon, const double Delta, const double x0, 
    # NumericVector theta, NumericMatrix eps)
    return(ts(rCIRESGcppexact(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                              theta = c(theta1, theta2, theta3), eps = eps), 
              start = 0, end = horizon, deltat = delta))
  }
  
  if (model == "OU")
  {
    if (is.null(theta1) || is.null(theta2) || is.null(theta3))
      stop("'theta1', 'theta2' and 'theta3' must be provided")
    
    if(length(theta1) != 1 || length(theta2) != 1 || length(theta3) != 1) 
      stop(paste("length required for 'theta1', 'theta2' and 'theta3' : 1"))
    
    # rOUESGcpp(const int N, const int horizon, const double Delta, const double x0, 
    # NumericVector theta, NumericMatrix eps)     
    return(ts(rOUESGcppexact(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                             theta = c(theta1, theta2, theta3), eps = eps), 
              start = 0, end = horizon, deltat = delta))
  }
  
  if (model == "GBM")
  {
    if (is.null(mu.z) && is.null(sigma.z) && is.null(p) && is.null(eta_up) && is.null(eta_down))
    {
      # rGBMESGcpp(const int N, const int horizon, const double Delta, const double x0, 
      # NumericMatrix theta1, NumericMatrix theta2, NumericMatrix eps)       
      return(ts(rGBMESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                           theta1 = theta1, theta2 = theta2, eps = eps), 
                start = 0, end = horizon, deltat = delta))
    }
    
    if (!is.null(lambda) && !is.null(mu.z) && !is.null(sigma.z))
    {
      # rGBMjumpsnormESGcpp(const int N, const int horizon, const double Delta,
      # const double x0, NumericMatrix theta1, NumericMatrix theta2,
      # const double lambda, const double mu, const double sigma, 
      # NumericMatrix eps)       
      return(ts(rGBMjumpsnormESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                                    theta1 = theta1, theta2 = theta2,
                                    lambda = lambda, mu = mu.z, sigma = sigma.z, 
                                    eps = eps), start = 0, end = horizon, deltat = delta))
    }
    
    if (!is.null(lambda) && !is.null(p) && !is.null(eta_up) && !is.null(eta_down))
    {
      # rGBMjumpskouESGcpp(const int N, const int horizon, const double Delta,
      # const double x0, NumericMatrix theta1, NumericMatrix theta2,
      # const double lambda, const double eta_up, const double eta_down,
      # const double p, NumericMatrix eps)
      return(ts(rGBMjumpskouESGcpp(N = n, horizon = horizon, Delta = delta, x0 = x0, 
                                   theta1 = theta1, theta2 = theta2,
                                   lambda = lambda, eta_up = eta_up, eta_down = eta_down,
                                   p = p, eps = eps), 
                start = 0, end = horizon, deltat = delta))
    }
  }
  }                                                                                                                                  
simdiff <- cmpfun(simdiff)


# Instantaneous forward rates ---------------------------------------------


#'@title 
#'
#'Instantaneous forward rates
#'
#'@description
#'
#'This function provides instantaneous forward rates. They can be used 
#'in no-arbitrage short rate models, to fit the yield curve exactly. 
#'
#'@param in.maturities input maturities
#'
#'@param in.zerorates input zero rates
#'
#'@param n number of independent observations
#'
#'@param horizon horizon of projection
#'
#'@param out.frequency either "annual", "semi-annual", "quarterly", "monthly", 
#'"weekly", or "daily" (1, 1/2, 1/4, 1/12, 1/52, 1/252)
#'
#'@param ... additional parameters provided to \code{\link{ycinter}}
#'
#'@author
#'
#'Thierry Moudiki
#'
#'@references
#'
#' Thierry Moudiki (2013). ycinterextra: Yield curve or zero-coupon prices interpolation and extrapolation. R package version 0.1. URL 
#'\url{http://CRAN.R-project.org/package=ycinterextra}
#'
#'@examples
#'
#'# Yield to maturities
#'txZC <- c(0.01422,0.01309,0.01380,0.01549,0.01747,0.01940,0.02104,0.02236,0.02348,
#'          0.02446,0.02535,0.02614,0.02679,0.02727,0.02760,0.02779,0.02787,0.02786,0.02776
#'          ,0.02762,0.02745,0.02727,0.02707,0.02686,0.02663,0.02640,0.02618,0.02597,0.02578,0.02563)
#'
#'# Observed maturities
#'u <- 1:30
#'
#'\dontrun{
#'par(mfrow=c(2,2))
#'fwdNS <- esgfwdrates(in.maturities = u, in.zerorates = txZC, 
#'                     n = 10, horizon = 20, 
#'                     out.frequency = "semi-annual", method = "NS")
#'matplot(time(fwdNS), fwdNS, type = 'l')
#'
#'fwdSV <- esgfwdrates(in.maturities = u, in.zerorates = txZC, 
#'                     n = 10, horizon = 20, 
#'                     out.frequency = "semi-annual", method = "SV")
#'matplot(time(fwdSV), fwdSV, type = 'l')
#'
#'fwdSW <- esgfwdrates(in.maturities = u, in.zerorates = txZC, 
#'                     n = 10, horizon = 20, 
#'                     out.frequency = "semi-annual", method = "SW")
#'matplot(time(fwdSW), fwdSW, type = 'l')
#'
#'fwdHCSPL <- esgfwdrates(in.maturities = u, in.zerorates = txZC, 
#'                        n = 10, horizon = 20, 
#'                        out.frequency = "semi-annual", method = "HCSPL")
#'matplot(time(fwdHCSPL), fwdHCSPL, type = 'l')
#'}
esgfwdrates <- function(in.maturities, in.zerorates,
                        n, horizon, 
                        out.frequency = c("annual", "semi-annual", 
                                          "quarterly", "monthly", 
                                          "weekly", "daily"),  
                        ...)
{
  if(is.null(in.maturities) || is.null(in.zerorates))
    stop("Zero rates and maturities must be provided")
  
  if (floor(n) != n) 
    stop("'n' must be an integer")
  
  if (floor(horizon) != horizon) 
    stop("'horizon' must be an integer")
  
  max.in.maturities <- max(in.maturities)
  
  if (horizon > max.in.maturities)
    stop("not enough maturities for the provided horizon")
  
  frequency <- match.arg(out.frequency)
  delta <- switch(frequency, 
                  "annual" = 1, 
                  "semi-annual" = 0.5, 
                  "quarterly" = 0.25, 
                  "monthly" = 1/12,
                  "weekly" = 1/52,
                  "daily" = 1/252)
  
  # simply coumpounded zero-coupon price
  pricefromeuribor <- function(t, T, L)
  {
    # Brigo P. 7
    tau <- T-t
    return(as.vector(1/(1 + L*tau)))    
  } 
  
  p <- pricefromeuribor(t = 0, T = in.maturities, L = in.zerorates)
  tt <- seq(from = min(in.maturities), to = max.in.maturities, by = delta)
  nn <- length(tt)
  
  yc <- ycinter(matsin = in.maturities, matsout = tt, p = p, 
                      typeres="prices", ...)  
  
  ZC.prices <- fitted(yc)
  
  ZC.prices <- c(1 + ((p[1] - 1)/(in.maturities[1] - 0))*(seq(0, 1, by = delta) - 0), 
                 ZC.prices[-1])
  
  tt <- seq(0, max(in.maturities), by = delta)
  
  nn <- length(tt)
  
   fwdrates <- ts(data = (1/delta)*(ZC.prices[-nn]/ZC.prices[-1]-1), 
                  start = 0,
                  deltat = delta)     
   
   return(ts(replicate(n, fwdrates), start = 0,
             deltat = delta))
}
esgfwdrates <- cmpfun(esgfwdrates)