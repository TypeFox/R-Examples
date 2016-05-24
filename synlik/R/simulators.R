############################################################################################################################################
##################################                    S I M U L A T O R S                       ############################################
############################################################################################################################################
# All simulators and wrappers around C/C++ simulators are here

########
# Blowfly full model (Simon's simulators) _INTERNAL_
########
## lu,lu1 are noise matrices of the same dimension.
## they contain ilogistic(uniform) random deviates,
## fixed from rep to rep to avoid monte carlo error. 
## ncol(lu) is number of reps, nrow(lu) is length of sim 
## including burn.in period, theta is parameter vector
## first transform the ilogistic(uniform) deviates in lu and lu1 to
## gamma deviates with mean 1 and given variance, and do it quickly!
.blowfly <- function (theta, lu, lu1, burn.in, pop.siz) 
{
  ## lu,lu1 are noise matrices of the same dimension.
  ## they contain logit(uniform) random deviates,
  ## fixed from rep to rep to avoid monte carlo error. 
  ## ncol(lu) is number of reps, nrow(lu) is length of sim 
  ## including burn.in period, theta is parameter vector
  
  ## first transform the logit(uniform) deviates in lu and lu1 to
  ## gamma deviates with mean 1 and given variance, and do it quickly!
  var <- theta[4];
  var1 <- theta[6]
  
  luk <- seq(min(lu),max(lu),length=200) ## even spacing on logit(prob) scale
  uk <- exp(luk);uk <- uk/(1+uk)          ## on prob scale
  gk <- qgamma(uk,shape=1/var,scale=var) ## reference quantiles
  gk[!is.finite(gk)] <- 0                ## can underflow badly at lower end
  e <- approx(x=luk,y=gk,xout=lu)$y      ## fast approximate version of qgamma(u,...) 
  
  luk <- seq(min(lu1),max(lu1),length=200)
  uk <- exp(luk);uk <- uk/(1+uk)
  gk <- qgamma(uk,shape=1/var1,scale=var1)
  gk[!is.finite(gk)] <- 0 
  e1 <- approx(x=luk,y=gk,xout=lu1)$y 
  
  n.t <- nrow(lu)-burn.in
  n.reps <- ncol(lu)
  n <- matrix(0,n.t,n.reps)
  oo <- .C("blowC", n=as.double(n), as.double(theta), as.double(e), as.double(e1), as.integer(burn.in), as.integer(n.t),
           as.integer(n.reps), PACKAGE="synlik")
  t(matrix(oo$n, n.t, n.reps))
}

################
########
#' Simulates from the ricker model
#' 
#' @description Simulator for the stochastic Ricker model, as described by Wood (2010). The observations are Y_t ~ Pois(Phi * N_t),
#'              and the dynamics of the hidden state are given by N_t = r * N_\{t-1\} * exp( -N_\{t-1\} + e_t ), where e_t ~ N(0, Sigma^2).
#' 
#' @param param vector of log-parameters: logR, logSigma, logPhi. Alternatively a matrix \code{nsim} by 3 were each row is
#'              a different parameter vector.
#' @param nsim Number of simulations from the model.
#' @param extraArgs A named list of additional arguments: \itemize{
#'  \item{\code{nObs} = Length of each simulated time series.}
#'  \item{\code{nBurn} = Number of initial steps to be discarded before saving the following \code{nObs} steps.}
#'  \item{\code{randInit} = if \code{TRUE} (default) the initial state N0 is \code{runif(0, 1)}, otherwise it is equal to \code{extraArgs$initVal}.}
#'  \item{\code{initVal} = initial value N0, used only if \code{extraArgs$randInit == TRUE}.}
#' } 
#'
#' @param ... Need for compatibility with \code{synlik}, but not used.
#'
#' @return A matrix \code{nsim} by \code{nObs}, where each row is a simulated path.
#' 
#' @references  Simon N Wood. Statistical inference for noisy nonlinear ecological dynamic systems. Nature, 466(7310):1102--1104, 2010. \cr \cr       
#' @author Simon Wood and Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @seealso \link{ricker_sl}
#' @examples
#' 
#' tmp <- rickerSimul(c(3.8, -1.2, 2.3), nsim = 2, extraArgs = list("nObs" = 50, "nBurn" = 200))
#' matplot(t(tmp), type = 'l', ylab = "Y", xlab = "Time")
#' 
#' parMat <- rbind(c(3.8, -1.2, 2.3),  # Chaotic
#'                 c(2.5, -1.2, 2.3))  # Not Chaotic
#'                 
#' par(mfrow = c(2, 1))
#' tmp <- rickerSimul(parMat, nsim = 2, extraArgs = list("nObs" = 50, "nBurn" = 200))
#' plot(tmp[1, ], type = 'l', ylab = "Y", xlab = "Time")
#' plot(tmp[2, ], type = 'l', ylab = "Y", xlab = "Time")
#' @export
rickerSimul <- function(param, nsim, extraArgs, ...)
{
  if( !all( c("nObs", "nBurn") %in% names(extraArgs) ) ) stop("extraArgs should contain nBurn and nObs")
  nBurn <- extraArgs$nBurn
  nObs <- extraArgs$nObs
  
  if( is.null(extraArgs$randInit) ) randInit = TRUE else randInit <- extraArgs$randInit
  
  if( is.null(extraArgs$initVal) ) initVal = 1.0 else initVal <- extraArgs$initVal
  
  if( !is.matrix(param) ) param <- matrix(param, 1, length(param))
  
  .Call( "simpleModelsWrap", model = "ricker", days = nObs, nSimul = nsim, param = param, nBurn = nBurn, randInit = randInit, initVal = initVal, PACKAGE = "synlik" )
}


########
#' Simulates from the blowfly model
#' 
#' @description Simulator for the blowfly model proposed by Wood (2010).
#' 
#' @param param vector of log-parameters: delta, P, N0, var.p, tau and var.d. The interpretation of these parameters is
#'              described in Wood (2010).
#' @param nsim Number of simulations from the model.
#' @param extraArgs A named list of additional arguments: \itemize{
#'  \item{\code{nObs} = Length of each simulated time series.}
#'  \item{\code{nBurn} = Number of initial steps to be discarded before saving the following \code{nObs} steps.}
#'  \item{\code{steps} = Positive integer. If \code{steps == n} the observations correspond to \code{n} time steps.} 
#' } 
#' @param ... Need for compatibility with \code{synlik}, but not used.
#'
#' @return A matrix \code{nsim} by \code{nObs}, where each row is a simulated path.
#' 
#' @references  Simon N Wood. Statistical inference for noisy nonlinear ecological dynamic systems. Nature, 466(7310):1102--1104, 2010. \cr \cr
#'              Brillinger, D. R., J. Guckenheimer, P. Guttorp, and G. Oster. 1980. 
#'              Empirical modelling of population time series data: the case of age and density dependent 
#'              vital rates. Lectures on Mathematics in the Life Sciences13:65-90.  \cr \cr
#'              Nicholson, A. J. 1957. The self-adjustment of populations to change. 
#'              Cold Spring Harbor Symposia on Quantitative Biology22:153-173.        
#' @author Simon Wood and Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @seealso \link{blow_sl}
#' @examples
#' tmp <- blowSimul(param = log( c( "delta" = 0.16, "P" = 6.5, "N0" = 400, 
#'                                  "var.p" = 0.1, "tau" = 14, "var.d" = 0.1)  ), 
#'                  nsim = 2, 
#'                  extraArgs = list("nObs" = 200, "nBurn" = 1000, "steps" = 2))
#' matplot(t(tmp), type = 'l', ylab = "Y", xlab = "Time")
#' @export

blowSimul <- function(param, nsim, extraArgs, ...)
{
  if( !is.vector(param) ) param <- as.vector(param)
  if( length(param) != 6) stop("Wrong number of parameters")
  
  if( !all( c("nObs", "nBurn", "steps") %in% names(extraArgs) ) ) stop("extraArgs should contain nBurn, nObs and steps")
  nBurn <- extraArgs$nBurn
  nObs <- extraArgs$nObs
  steps <- extraArgs$steps
  
  totStep <- nBurn + nObs * steps
  
  noise1 <- matrix(.ilogistic(runif(nsim * totStep)), totStep, nsim)
  noise2 <- matrix(.ilogistic(runif(nsim * totStep)), totStep, nsim)
  
  simul <- .blowfly(theta = exp(param), 
                    lu = noise1, 
                    lu1 = noise2, 
                    burn.in = nBurn, pop.siz = nsim)[ , cumsum(c(1,rep(steps, nObs-1)))]
  
  return( simul )
}














