#######################
###### Function that checks the multivariate normality of the statistics
#######################
#' Checking the multivariate normal approximation.
#' @description Given an object of class \code{synlik} this routine provides a 
#'              graphical check of whether the distribution of the random 
#'              summary statistics is multivariate normal.
#' @param object An object of class \code{synlik} or a matrix where each row is a random vector.
#' @param param A vector of model's parameters at which the summary statistics will be simulated.
#' @param observed A vector of observed summary statistics. By default \code{NULL}, so \code{object@@data} will be used as observed statistics.
#'                 It will be looked at only if \code{object} is a matrix.
#' @param nsim number of summary statistics to be simulated if object is of class \code{synlik}, otherwise
#'             it is not used.
#' @param cex.axis Axis scale expansion factor.
#' @param cex.lab  Axis label expansion factor.
#' @param ... additional arguments to be passed to \code{object@@simulator} and \code{object@@summaries}.
#'            In general I would avoid using it and including in \code{object@@extraArgs} everything they need.
#' @details The method is from section 7.5 of Krzanowski (1988). The replicate vectors of summary statistic \code{S} 
#'          are transformed to variables which should be univariate chi squared r.v.s with DoF given by the number of rows of \code{S}. 
#'          An appropriate QQ-plot is produced, and the proportion of the data differing substantially from the ideal 
#'          line is reported. Deviations at the right hand end of the plot indicate that the tail behaviour of the Normal 
#'          approximation is poor: in the context of synthetic likelihood this is of little consequence.  
#'          Secondly, \code{s} is transformed to a vector which should be i.i.d. N(0,1) under
#'          multivariate normality, and a QQ plot is produced. Unfortunately this approach is not very
#'           useful unless the dimension of \code{s} is rather large. In simulations, perfectly MVN data produce 
#'           highly variable results, so that the approach lacks any real power. 
#' @return  Mainly produces plots and prints output. Also an array indicating 
#'          proportion of simulated statistics smaller than observed.
#' @references Krzanowski, W.J. (1988) Principles of Multivariate Analysis. Oxford.            
#' @author Simon N. Wood, maintained by Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @examples
#' #### Create Object
#' data(ricker_sl)
#' 
#' #### Simulate from the object
#' ricker_sl@@data <- simulate(ricker_sl)
#' ricker_sl@@extraArgs$obsData <- ricker_sl@@data 
#' 
#' #### Checking multivariate normality
#' checkNorm(ricker_sl)
#' 
#' # With matrix input
#' checkNorm(matrix(rnorm(200), 100, 2))
#' @export
#'


checkNorm <- function(object,
                      param = object@param, 
                      nsim = 1e3,
                      observed = NULL,
                      cex.axis = 1, 
                      cex.lab = 1, 
                      ...)
{
  oldPar <- par(no.readonly = TRUE)
  
  if( is(object, "synlik") )
  {
  
  S <- t( simulate(object, nsim, param, stats = TRUE, ...) )
  
  # I copy this function so I can mtrace() it
  summaries <- object@summaries
  
  if( !is.null(object@data) )
  { 
    
    s <- if( !is.null(summaries) )
    {
      drop( summaries(x = object@data, extraArgs = object@extraArgs, ...) )
    } else
    {
      drop( object@data )
    }
    
  }else
  { 
    s <- NULL
  }
  
  } else{
    if( !is.numeric(object) ) stop("object should be either of class \"synlik\", a numeric matrix or vector")
    if( !is.matrix(object) ) object <- matrix(object, length(object), 1)
    S <- t( object )
    s <- drop( observed )
  }
  
  p <- nrow(S)
  n <- ncol(S)
  if (n < 10 * p) 
    warning("You don't really have enough repetitions (nsim) for this approach")
  ps <- s
  for (i in 1:nrow(S)) ps[i] <- sum(S[i, ] < s[i])/ncol(S)
  um <- robCov(S)
  ms <- as.numeric(um$mY)
  S <- S - ms
  z <- colSums(S * (t(um$E) %*% um$E %*% S))
  q <- log(qchisq((1:n - 0.5)/n, df = p))
  z <- log(sort(z))
  
  par(mfrow = c(2, 2))
  plot(q, z, type = "l", col = "grey", xlab = "log theoretical quantiles", 
       ylab = "log observed quantiles", cex.axis = cex.axis, 
       cex.lab = cex.lab)
  points(q, z, pch = ".")
  abline(0, 1, col = 2)
  cat("\nproportion |log(z)-log(q)|>.25 = ", sum(abs(z - q) > 
                                                   0.25)/n, "\n")
  if (!is.null(s)) {
    z <- um$E %*% (s - ms)
    q <- sum(z^2)
    abline(h = log(q), lty = 2)
  }
  for (i in 1:nrow(S)) S[i, ] <- S[i, ]/um$sd[i]
  n <- ncol(S)
  z <- qnorm((1:n - 0.5)/n)
  rz <- range(z)
  rz <- c(rz[1] - 2, rz[2] + 2)
  plot(z, sort(S[1, ]), type = "l", col = "grey", ylim = rz, 
       xlab = "N(0,1) quantiles", ylab = "marginal quantiles", 
       main = "Marginal Q-Q Plot",
       cex.axis = cex.axis, cex.lab = cex.lab)
  points(z, sort(S[1, ]), pch = ".")
  
  if(nrow(S) > 1)
  {
    for (i in 2:nrow(S)) 
    {
      lines(z, sort(S[i, ]), col = "grey")
      points(z, sort(S[i, ]), pch = ".")
    }
  }
  abline(0, 1, col = 2)
  
  if (!is.null(s)) {
    z <- um$E %*% (s - ms)
    qqnorm(z, cex.axis = cex.axis, cex.lab = cex.lab)
    qqline(z, col = 2)
  }
  
  par(oldPar)
  
  return(ps)
}
