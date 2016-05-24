#' Calculate Time-dependent ROC and AUC
#'
#' This is the main function of this package. It calculates the time-dependent
#'  sensitivity and specificity and area under the curve (AUC) using a nonparametric weighting adjustment.
#'  It also provides variance estimation through bootstrap.
#'
#'@param X a numeric vector of biomarker values. Same length with \code{Y} and \code{delta}.
#'@param Y a numeric vector of time to event.Same length with \code{X} and \code{delta}.
#'@param delta a vector of binary indicator of event (1) or censoring (0). Same length with \code{X} and \code{Y}.
#'@param tau a scalar, the prediction horizon at which the prediction is evaluated.
#'@param span a numeric value, the proportion of neighbour observations used in nearest neighbor method, default is 0.1.
#'@param h a numeric value, the bandwidth of kernel weights, defualt is \code{NULL}. If not specified, the function will use the value of
#'          \code{span} to calculate kernel weights. In case both \code{span} and \code{h} are specified, the function will use \code{h}.
#'@param type a character value, indicating the type of kernel function used to calculate kernel weights. Default is "\code{uniform}" kernel. Other options are "\code{Epanechnikov}" and "\code{normal}".
#'        It will only be used when the bandwidth \code{h} is specified.
#'@param nboot the number of bootstrap replications to be used for variance estimation; default is \code{nboot = 0}, corresponding to no variance estimation.
#'@param alpha \code{1-level of confidence interval}, default is \code{0.05}. It is used only when \code{nboot > 0}.
#'@param n.grid an positive integer, the number of grid points used when calculating the ROC curve. The default is \code{1000}.
#'@param cut.off a vector of biomarker cut-off values at which sensitivity and specificity will be calculated.When bootstrap is requested, the corresponding confidence intervals will also be provided.
#'@param X.min the lower boundary of grid cut-off points for biomarker \code{X}. If not specified, default will be the minimum of biomarker \code{X}.
#'@param X.max the upper boundary of grid cut-off points for biomarker \code{X}. If not specified, default will be the maximum of biomarker \code{X}.
#'        \code{X.min} and \code{X.max} are not needed for point estimate but are needed for bootstrap variance estimation.
#'@details This function read in the prognostic biomarker value \code{X}, the time-to-event data \code{Y} and censoring indicator \code{delta} to calculate
#'  the sensitivity and specificity at the prediction horizon \code{tau} for a series specified grid points. It uses a simple
#'  nonparametric weight adjustments for right censored data (Li \emph{et al.}, 2015).
#'@seealso \code{\link[survival]{survfit}},  \code{\link[survivalROC]{mayo}}
#'@return Returns a list of the following items:
#'
#'@return \code{ROC:} a data frame of dimension \code{(2+n.grid) x 3}, the three columns are: \code{grid}, \code{sens}, and \code{spec}.
#'@return \code{AUC:} a data frame of one row and four columns: \code{AUC}, standard error of \code{AUC}, the lower and upper limits of bootstrap CI.
#'  \code{AUC} is calculated by integrating the area under ROC curve with trapezoidal method.
#'@return \code{AUC2:} a data frame of one row and four columns: \code{AUC2}, standard error of \code{AUC2}, the lower and upper limits of bootstrap CI.
#'  \code{AUC2} is the AUC calculated by the concordance based formula (Li \emph{et al.}, 2015).
#'@return \code{prob:} a data frame of three columns if \code{nboot=0}: \code{cut.off}, \code{sens}, and \code{spec}. If \code{nboot>0}, another six
#'        columns of standard error, lower and upper limits of both \code{sens} and \code{spec} will be added. The number of rows equals length of \code{cut.off}.
#'              A series of sensivitity and specificity are calculated at requested \code{cut.off} points.
#'@importFrom survival survfit
#'@examples
#'library( survival ) ;
#'data( mayo ) ;
#'dat <- mayo[ ,c( "time","censor","mayoscore5" )] ;
#'
#'fm <- tdROC( X = dat$mayoscore5, Y = dat$time, delta = dat$censor,
#'        tau = 365*6, span = 0.1, nboot = 0, alpha = 0.05,
#'        n.grid = 1000, cut.off = 5:9 ) ;
#'
#'@export
#'@author Liang Li, Cai Wu
#'@references Li, Liang, Bo Hu, and Tom Greene. "A Simple Method to Estimate the Time-dependent ROC Curve Under Right Censoring." (2015).
#'            \url{http://biostats.bepress.com/cobra/art114/}


tdROC <- function( X, Y, delta, tau, span = 0.1, h=NULL, type="uniform", cut.off = NULL,
                   nboot = 0, alpha = 0.05, n.grid = 1000,
                   X.min = NULL, X.max = NULL ) {
  # Calculate the time-dependent sensitivity and specificity and
  # area under the curve (AUC)
  # Arguments:
  #  -- X: the vector of biomarker values
  #  -- Y: time to event
  #  -- delta: indicator of event (1) or censoring (0)
  #  -- tau: the prediction horizon
  #  -- span: the proportion of observations used in calculating kernel weights
  #           (i.e., bandwidth in nearest neighbor method)
  #  -- nboot: number of bootstrap, to be used the variance estimation;
  #            nboot = 0 corresponds to no variance estimation
  #  -- alpha: 1-level of confidence interval, default to be 0.05. It is
  #            used only when nboot > 0
  #  -- n.grid: number of biomarker cut off values used
  #             when calculating the ROC curve
  #  -- cut.off: a vector of biomarker cut.off values at which sensitivity and
  #              specificity will be calculated
  #  -- X.min, X.max: the min and max of X; if \code{NULL}, they will be calculated
  #             inside the function; these are not needed for point estimate,
  #             they are needed for bootstrap
  # Return:
  #  -- ROC: a data frame of three columns:
  #             grid, sensitivity and specificity
  #  -- AUC: a data frame of one row and four columns
  #             AUC, standard error of AUC, the lower and upper bootstrap CI
  #             AUC2 is the AUC calculated by integrating area under ROC curve
  #  -- AUC2: a data frame of one row and four columns
  #             AUC2, standard error of AUC2, the lower and upper bootstrap CI
  #             AUC2 is the AUC calculated by the concordance based formula
  #  -- prob: a data frame of three columns:
  #             cut.off, sensitivity and specificity
  # NOTE: X, Y and delta must have the same length
  n <- length(X) ;
  positive <- rep(NA, n) ;
  for (i in 1:n) {
    if ( Y[i] > tau ) {
      positive[i] <- 0 ;
    } else {
      if ( delta[i] == 1 ) {
        positive[i] <- 1 ;
      } else {
        kw <- calc.kw( X=X, x0=X[i], span=span, h=h, type=type ) ;
        fm <- survfit( Surv(Y, delta) ~ 1, weights = kw ) ;
        tmp <- summary(fm, times = c(Y[i], tau))$surv ;
        if ( tmp[1] == 0 ) {
          positive[i] <- 1 ;
        } else {
          positive[i] <- 1 - tmp[2]/tmp[1] ;
        }
      }
    }
  }
  negative <- 1 - positive ;
  if ( is.null(X.min) ) { X.min <- min(X) }
  if ( is.null(X.max) ) { X.max <- max(X) }
  grid <- c( -Inf, seq( X.min, X.max, length=n.grid ), Inf ) ;
  sens <- spec <- NULL ;
  for (this.c in grid ) {
    sens <- c( sens, sum(positive*as.numeric(X > this.c))/sum(positive) ) ;
    # sensitivity that incorporates fractional "positive"
    spec <- c( spec, sum(negative*as.numeric(X <= this.c))/sum(negative) ) ;
    # specificity that incorporates fractional "negative"
  }
  ROC <- data.frame( grid = grid,
                     sens = sens,
                     spec = spec ) ;
  # calculate the first estimator for AUC using integral
  AUC <- data.frame( value = calc.AUC( sens, spec ),
                     sd = NA,
                     lower = NA,
                     upper = NA ) ;
  # calculate the second estimator for AUC using formula
  W <- positive ;
  numer <- denom <- 0 ;
  for (i in 1:n) {
    for (j in 1:n) {
      numer <- numer + W[i]*(1-W[j])*( as.numeric(X[i] > X[j]) +
                                         0.5*as.numeric(X[i] == X[j]) ) ;
      denom <- denom + W[i]*(1-W[j]) ;
    }
  }
  AUC2 <- data.frame( value = numer/denom ,
                      sd = NA, lower = NA, upper = NA ) ;
  # calculate the sensitivity and specificity at selected
  # (specified in the arguments) cut-offs
  sens <- spec <- NULL ;
  if ( !is.null(cut.off) ) {
    for (this.c in cut.off ) {
      sens <- c( sens, sum(positive*as.numeric(X > this.c))/sum(positive) ) ;
      # sensitivity that incorporates fractional "positive"
      spec <- c( spec, sum(negative*as.numeric(X <= this.c))/sum(negative) ) ;
      # specificity that incorporates fractional "negative"
    }
    prob <- data.frame( cut.off = cut.off,
                        sens = sens,
                        spec = spec ) ;
  } else {
    prob <- NULL ;
  }

  if ( nboot > 0 ) {
    # start bootstrap for AUC
    boot.AUC <- boot.AUC2 <- rep(NA, nboot) ;
    if ( !is.null(cut.off) ) {
      boot.sens <- matrix( NA, nrow=nboot, ncol=length(cut.off) ) ;
      boot.spec <- matrix( NA, nrow=nboot, ncol=length(cut.off) ) ;
    }
    set.seed(123) ;
    # the random number seed is hardcoded
    for (b in 1:nboot) {
      loc <- sample( x = 1:n, size = n, replace = T ) ;
      X2 <- X[loc] ;
      Y2 <- Y[loc] ;
      delta2 <- delta[loc] ;
      out <- tdROC( X2, Y2, delta2, tau, span, nboot = 0, alpha, n.grid,
                    cut.off = cut.off, X.min = X.min, X.max = X.max ) ;
      boot.AUC[b] <- out$AUC$value ;
      boot.AUC2[b] <- out$AUC2$value ;
      if ( !is.null(cut.off) ) {
        boot.sens[b, ] <- out$prob$sens ;
        boot.spec[b, ] <- out$prob$spec ;
      }
    }
    tmp1 <- sd(boot.AUC) ;
    tmp2 <- as.numeric( quantile( boot.AUC, prob = c(alpha/2, 1-alpha/2) ) ) ;
    AUC$sd <- tmp1 ;
    AUC$lower <- tmp2[1] ;
    AUC$upper <- tmp2[2] ;
    #
    tmp1 <- sd(boot.AUC2) ;
    tmp2 <- as.numeric( quantile( boot.AUC2, prob = c(alpha/2, 1-alpha/2) ) ) ;
    AUC2$sd <- tmp1 ;
    AUC2$lower <- tmp2[1] ;
    AUC2$upper <- tmp2[2] ;
    #
    if ( !is.null(cut.off) ) {
      prob$sens.sd <- apply( boot.sens, 2, sd ) ;
      prob$sens.lower <- apply( boot.sens, 2, quantile, prob = alpha/2 ) ;
      prob$sens.upper <- apply( boot.sens, 2, quantile, prob = 1-alpha/2 ) ;
      prob$spec.sd <- apply( boot.spec, 2, sd ) ;
      prob$spec.lower <- apply( boot.spec, 2, quantile, prob = alpha/2 ) ;
      prob$spec.upper <- apply( boot.spec, 2, quantile, prob = 1-alpha/2 ) ;
    } else {
      prob$sens.sd <- NA ;
      prob$sens.lower <- NA ;
      prob$sens.upper <- NA ;
      prob$spec.sd <- NA ;
      prob$spec.lower <- NA ;
      prob$spec.upper <- NA ;
    }
  }

  pct.ctrl <- mean( Y > tau ) ;
  pct.case <- mean( Y <= tau & delta == 1 ) ;
  pct.not.sure <- mean( Y <= tau & delta == 0 ) ;
  return( list( ROC = ROC, AUC = AUC, AUC2 = AUC2, prob = prob ) ) ;
}
