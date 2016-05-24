#' Simulates statistics related to exceedance region.
#' 
#' \code{statistic.sim} simulates statistics related to the construction of confidence regions for exceedance sets and contour lines.
#' 
#' When \code{alternative = "two.sided"}, the \code{...} argument must include \code{user.cov} (a user-specified covariance function), \code{pgrid} (the grid of locations to be predicted, produced by \code{create.pgrid} or \code{create.pgrid2}), \code{X} (the matrix of covariates for the observed data), and any other arguments needed by \code{user.cov}. Note that \code{user.cov} should take \code{cLcoords} as its first argument (a matrix containing the coordinates of contour lines under consideration). Additional arguments to \code{user.cov} are passed internally using the \code{...} argument. The \code{user.cov} function should return a list with values \code{V} (the covariance matrix of the observed data), \code{Vop} (the cross-covariance matrix between the observed data and the responses with coordinates in cL), \code{Vp} (the covariance matrix of the responses with coordinates in \code{cL}), and \code{Xp} (the matrix of covariates for the coordinates contained in \code{cL}). See the Examples section.
#' 
#' @param krige.obj An object from the function \code{krige.uk} in the \code{SpatialTools} package.
#' @param level The threshold/exceedance level under consideration.
#' @param alternative Indicates the type of exceedance region or level curve under consideration.  For exceedances above a threshold, use (\code{alternative = "less"}).  For exceedances below a threshold, use (\code{alternative = "greater"}).  For contour lines, use (\code{alternative = "two.sided"}). Defaults to "less".
#' @param ... Additional arguments when \code{alternative = "two.sided"}. See Details.
#' 
#' @return Returns a list with components: 
#' \item{statistic}{A vector with the observed values of the test statistic.}
#' \item{statistic.sim}{A vector with the observed values of the test statistic.}
#' \item{alternative}{The alternative hypothesis provided to \code{statistic.sim}.}
#' \item{level}{The threshold level under consideration.}
#' 
#' @author Joshua French
#' @import SpatialTools
#' @export
#' @examples
#' library(SpatialTools)
#' 
#' # Example for exceedance regions
#' 
#' set.seed(10)
#' # Load data
#' data(sdata)
#' # Create prediction grid
#' pgrid <- create.pgrid(0, 1, 0, 1, nx = 26, ny = 26)
#' pcoords <- pgrid$pgrid
#' # Create design matrices
#' coords = cbind(sdata$x1, sdata$x2)
#' X <- cbind(1, coords)
#' Xp <- cbind(1, pcoords)
#' 
#' # Generate covariance matrices V, Vp, Vop using appropriate parameters for
#' # observed data and responses to be predicted
#' spcov <- cov.sp(coords = coords, sp.type = "exponential", sp.par = c(1, 1.5),
#'  error.var = 1/3, finescale.var = 0, pcoords = pcoords)
#' 
#' # Predict responses at pgrid locations
#' krige.obj <- krige.uk(y = as.vector(sdata$y), V = spcov$V, Vp = spcov$Vp, 
#'  Vop = spcov$Vop, X = X, Xp = Xp, nsim = 50, 
#'  Ve.diag = rep(1/3, length(sdata$y)) , method = "chol")
#'                 
#' # Simulate distribution of test statistic for different alternatives
#' statistic.sim.obj.less <- statistic.sim(krige.obj = krige.obj, level = 5,
#'  alternative = "less")
#' statistic.sim.obj.greater <- statistic.sim(krige.obj = krige.obj, level = 5,
#'  alternative = "greater")
#' # Construct null and rejection sets for two scenarios
#' n90 <- exceedance.ci(statistic.sim.obj.less, conf.level = .90, type = "null")
#' r90 <- exceedance.ci(statistic.sim.obj.greater,conf.level = .90, 
#'  type = "rejection")       
#' # Plot results
#' plot(pgrid, n90, col="blue", add = FALSE, xlab = "x", ylab = "y")
#' plot(pgrid, r90, col="orange", add = TRUE)
#' legend("bottomleft", 
#'  legend = c("contains true exceedance region with 90 percent confidence", 
#'  "is contained in true exceedance region with 90 percent confidence"),
#'  col = c("blue", "orange"), lwd = 10)  
#' 
#' # Example for level curves
#' data(colorado)
#' ocoords <- colorado$ocoords
#' odata <- colorado$odata    
#' 
#' # Set up example
#' nsim <- 50
#' u <- log(16)
#' np <- 26
#' conf.level <- 0.90   
#' x.min <- min(ocoords[,1])
#' x.max <- max(ocoords[,1])
#' y.min <- min(ocoords[,2])
#' y.max <- max(ocoords[,2])   
#' 
#' #pixelize the domain
#' pgrid <- create.pgrid(x.min, x.max, y.min, y.max, nx = np, ny = np)
#' pcoords <- pgrid$pgrid; upx <- pgrid$upx; upy <- pgrid$upy
#' names(pcoords) <- c("lon", "lat")    
#' 
#' # Set up covariates matrices
#' X <- cbind(1, ocoords)
#' Xp <- cbind(1, pcoords)           
#' 
#' # Estimate covariance parameters
#' cov.est <- maxlik.cov.sp(X, odata, sp.type = "exponential", range.par = 1.12,
#'  error.ratio = 0.01, reml = TRUE, coords = ocoords)
#' 
#' # Create covariance matrices
#' myCov <- cov.sp(coords = ocoords, sp.type = "exponential", 
#'  sp.par = cov.est$sp.par, error.var = cov.est$error.var, pcoords = pcoords)
#' 
#' # Krige and do conditional simulation 
#' krige.obj <- krige.uk(y = odata, V = myCov$V, Vp = myCov$Vp, Vop = myCov$Vop, 
#'  X = X, Xp = Xp, nsim = nsim, Ve.diag = rep(cov.est$error.var, 
#'  length(odata)))
#' 
#' # Create user covariance function for simulating statistic for confidence 
#' # regions
#' user.cov <- function(cLcoords,...)
#' {
#'    arglist <- list(...)
#'    coords <- arglist$coords
#'    sp.type <- arglist$sp.type
#'    sp.par <- arglist$sp.par
#'    V <- arglist$V
#'    out <- list(V = arglist$V,
#'                Vp = sp.par[1] * exp(-dist1(cLcoords)/sp.par[2]),
#'                Vop = sp.par[1] * exp(-dist2(coords, cLcoords)/sp.par[2]))
#'    out$Xp <- cbind(1, cLcoords)
#'    return(out)
#' }
#' 
#' # Simulation statistic for confidence regions
#' statistic.sim.obj <- statistic.sim(krige.obj = krige.obj, level = u, 
#'  alternative = "two.sided", user.cov = user.cov, y = odata, pgrid = pgrid, 
#'  X = X, coords = ocoords, pcoords = pcoords, V = myCov$V, 
#'  sp.type = "exponential", sp.par = cov.est$sp.par)
#' 
#' # Create 90% confidence region
#' n90 <- exceedance.ci(statistic.sim.obj, conf.level = conf.level, 
#'  type = "null")
#' # Get estimated contour lines
#' cL <- contourLines(pgrid$upx, pgrid$upy, matrix(krige.obj$pred, nrow = np), 
#'  level = u)
#' 
#' # Plot results
#' plot(ocoords, xlab = "longitude", ylab = "latitude", type = "n", 
#'  cex.lab = 1.5, cex.axis = 1.5)
#' plot(pgrid, n90, col = "grey", add = TRUE)
#' plot.contourLines(cL, col="black", lwd=2, lty = 2, add = TRUE)

statistic.sim <- function(krige.obj, level, alternative = "less", ...)
{
  if(is.null(krige.obj$sim))
  {
    stop("krige.obj$sim cannot be NULL.  Try setting nsim > 0.")
  }
  if(!is.numeric(level) || length(level) > 1)
  {
    stop("level must be a numeric vector of length 1")
  }
  if(!(alternative == "less" || alternative == "greater" || alternative == "two.sided"))
  {
    stop('alternative must equal "two.sided" or "less" or "greater"')
  }
  
  nsim <- ncol(krige.obj$sim)
  stat.sim <- numeric(nsim)
  if(alternative != "two.sided")
  {
    statistic <- (krige.obj$pred - level)/sqrt(krige.obj$mspe)
  }else
  {
    statistic <- abs(krige.obj$pred - level)/sqrt(krige.obj$mspe)
  }
  
  if(alternative == "less")
  {
    for(i in 1:nsim)
    {
      which.exceedance.sim <- which(krige.obj$sim[, i] >= level)
      if(length(which.exceedance.sim) > 0)
      {
        stat.sim[i] <- min(statistic[which.exceedance.sim])
      }
    }
  }else if(alternative == "greater")
  {
    for(i in 1:nsim)
    {
      which.exceedance.sim <- which(krige.obj$sim[, i] <= level)
      if(length(which.exceedance.sim) > 0)
      {
        stat.sim[i] <- max(statistic[which.exceedance.sim])
      }
    }
  }
  else
  {
    # Check arguments.  Create unspecified arguments if needed.
    arglist <- list(...)
    argnames <- names(arglist)
    user.cov <- arglist$user.cov
    pgrid <- arglist$pgrid
    npx <- length(pgrid$upx)
    npy <- length(pgrid$upy)
    
    for(m in 1:nsim)
    {
      simmat <- matrix(krige.obj$sim[, m], nrow = npx, ncol = npy)
      cL <- contourLines(pgrid$upx, pgrid$upy, simmat, levels = level)
      
      if(length(cL) > 0)
      {
        obj <- user.cov(cLcoords = get.contours(cL), ...)
        krige.cL <- krige.uk(arglist$y, obj$V, obj$Vp, obj$Vop, arglist$X, obj$Xp)
        stat.sim[m] <- max(abs(krige.cL$pred - level)/sqrt(krige.cL$mspe))
      }
    }
  }
  
  return(list(statistic = statistic, statistic.sim = stat.sim,  
              alternative = alternative, level = level))
}
