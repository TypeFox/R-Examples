#---------------------------------------------------------------------
# l2boost (DISCRETE/HYBRID/FRIEDMAN/LARS)
# 
# notes:
# nu is ignored in LARS/LIMIT
# rescales x
# y is centered --- > mean is returned as part of the object
# !!! remove x-columns with NA's !!! CAUTION
#---------------------------------------------------------------------
#' @export l2boost
#' @name l2boost
#' 
l2boost <- function(x, ...)UseMethod("l2boost")
#
#' @title Generic gradient descent boosting method for linear regression.
#' 
#' @description Efficient implementation of Friedman's boosting algorithm  [Friedman (2001)] with L2-loss function and coordinate
#'  direction (design matrix columns) basis functions. This includes the elasticNet data augmentation of Ehrlinger and Ishwaran (2012), 
#'  which adds an L2-penalization (lambda) similar to the elastic net [Zou and Hastie (2005)].
#' 
#' @details
#'  The \code{\link{l2boost}} function is an efficient implementation of a generic boosting method [Friedman (2001)] for
#' linear regression using an L2-loss function. The basis functions are the column vectors of the design matrix. 
#' \code{\link{l2boost}} scales the design matrix such that the coordinate columns of the design correspond to the
#' gradient directions for each covariate. The boosting coefficients are equivalent to the gradient-correlation of each 
#' covariate. Friedman's gradient descent boosting algorithm proceeds at each step along the covariate direction closest
#' (in L2 distance) to the maximal gradient descent direction.
#' 
#' We include a series of algorithms to solve the boosting optimization. These are selected through the \emph{type} argument
#' \itemize{
#'  \item \emph{friedman} - The original, bare-bones l2boost (Friedman (2001)). This method takes a fixed step size of length
#'    \emph{nu}.
#' \item \emph{lars} - The l2boost-lars-limit (See Efron et.al (2004)). This algorithm takes a single step of the 
#'   optimal length to the critical point required for a new coordinate direction to become favorable. Although optimal
#'   in the number of steps required to reach the OLS solution, this method may be computationaly expensive for large p
#'   problems, as the method requires a matrix inversion to calculate the step length. 
#' \item \emph{discrete} - Optimized Friedman algorithm to reduce number of evaluations required 
#'   [Ehrlinger and Ishwaran 2012]. The algorithm dynamically determines the number of steps of length \emph{nu} to take along
#'   a descent direction. The discrete method allows the algorithm to take step sizes of multiples of \emph{nu} at any evaluation.
#' \item \emph{hybrid} - Similar to discrete, however only allows combining steps along the first descent direction. 
#'   \emph{hybrid} Works best if \emph{nu} is moderate, but not too small. In this case, Friedman's algorithm would take 
#'   many steps along the first coordinate direction, and then cycle when multiple coordinates have similar gradient 
#'   directions (by the L2 measure).
#' }
#' 
#' \code{\link{l2boost}} keeps track of all gradient-coorelation coefficients (\emph{rho}) at each iteration in addition to the maximal
#' descent direction taken by the method. Visuallizing these coefficients can be informative of the inner workings of gradient boosting 
#' (see the examples in the \code{\link{plot.l2boost}} method).
#' 
#' The \code{\link{l2boost}} function uses an arbitrary L1-regularization parameter (nu), and includes the elementary 
#' data augmentation of Ehrlinger and Ishwaran (2012), to add an L2-penalization (lambda) similar to the elastic net 
#' [Zou and Hastie (2005)]. The L2-regularization reverses repressibility, a condition where one variable acts as 
#' a boosting surrogate for other, possibly informative, variables. Along with the decorrelation 
#' effect, this \emph{elasticBoost} regularization circumvents L2Boost deficiencies in correlated settings. 
#' 
#' We include a series of S3 functions for working  with \code{\link{l2boost}} objects:
#' \itemize{
#' \item \code{\link{print}} (\code{\link{print.l2boost}}) prints a summary of the \code{\link{l2boost}} fit.
#' \item \code{\link{coef}} (\code{\link{coef.l2boost}}) returns the \code{\link{l2boost}} model regression coefficients at any point 
#' along the solution path. 
#' \item \code{\link{fitted}} (\code{\link{fitted.l2boost}}) returns the fitted \code{\link{l2boost}} response estimates (from
#' the training dataset) along the solution path. 
#' \item \code{\link{residuals}} (\code{\link{residuals.l2boost}}) returns the training set \code{\link{l2boost}} residuals along the
#' solution path.
#' \item \code{\link{plot}} (\code{\link{plot.l2boost}}) for graphing model cofficients of an \code{\link{l2boost}} object.
#' \item \code{\link{predict}} (\code{\link{predict.l2boost}}) for generating \code{\link{l2boost}} prediction estimates on possibly 
#' new test set observations.
#' }
#' A cross-validation method (\code{\link{cv.l2boost}}) is also included for L2boost and elasticBoost, for cross-validated error estimats 
#' and regularization parameter optimizations.
#'
#' @references Friedman J. (2001) Greedy function approximation: A gradient boosting machine. \emph{Ann. Statist.}, 29:1189-1232
#' @references Ehrlinger J., and Ishwaran H. (2012). "Characterizing l2boosting" \emph{Ann. Statist.}, 40 (2), 1074-1101
#' @references Zou H. and Hastie T (2005) "Regularization and variable selection via the elastic net"  \emph{J. R. Statist. Soc. B}, 67, Part 2, pp. 301-320
#' @references Efron B., Hastie T., Johnstone I., and Tibshirani R. (2004). "Least Angle Regression" \emph{Ann. Statist.} 32:407-499
#'
#' @usage \method{l2boost}{default}(x, y, M, nu, lambda, trace, type, qr.tolerance, eps.tolerance, ...)
#'
#' @param x design matrix of dimension n x p
#' @param y response variable of length n
#' @param formula an object of class \code{\link{formula}} 
#'     (or one that can be coerced to that class): a symbolic 
#'     description of the model to be fitted. The details of 
#'     model specification are given under \code{\link{formula}}.
#' @param data an optional data frame, list or environment 
#'    (or object coercible by \code{\link{as.data.frame}} to 
#'    a data frame) containing the variables in the model used in the
#'    \code{\link{formula}}.
#' @param M number of steps to run boost algorithm (M >1)
#' @param nu L1 shrinkage parameter (0 < nu <= 1)
#' @param lambda L2 shrinkage parameter used for elastic net boosting (lambda > 0 || lambda = NULL)
#' @param type Choice of l2boost algorithm from "discrete", "hybrid", "friedman","lars". See details below. (default "discrete")
#' @param qr.tolerance tolerance limit for use in \code{\link{qr.solve}} (default: 1e-30)
#' @param eps.tolerance dynamic step size lower limit (default: .Machine$double.eps)
#' @param trace show runtime messages (default: FALSE)
#' @param ... other arguments (currently unused)
#'
#' @return A "l2boost" object is returned, for which print, plot, predict, and coef methods exist.
#' \item{call}{the matched call.}
#' \item{type}{Choice of l2boost algorithm from "friedman", "discrete", "hybrid", "lars"}
#' \item{nu}{The L1 boosting shrinkage parameter value}    
#' \item{lambda}{The L2 elasticNet shrinkage parameter value}  
#' \item{x}{The training dataset}
#' \item{x.na}{Columns of original design matrix with values na, these have been removed from x}
#' \item{x.attr}{scale attributes of design matrix}
#' \item{names}{Column names of design matrix}
#' \item{y}{training response vector associated with x, centered about the mean value ybar}
#' \item{ybar}{mean value of training response vector}
#' \item{mjk}{measure to favorability. This is a matrix of size p by m. Each coordinate j has a measure at each step m}
#' \item{stepSize}{vector of step lengths taken (\code{NULL} unless \code{type = "lars"})}       
#' \item{l.crit}{vector of column index of critical direction} 
#' \item{L.crit}{number of steps along each l.crit direction}      
#' \item{S.crit}{The critical step value where a direction change occurs}
#' \item{path.Fm}{estimates of response at each step m}
#' \item{Fm}{estimate of response at final step M}
#' \item{rhom.path}{boosting parameter estimate at each step m}   
#' \item{betam.path}{beta parameter estimates at each step m. List of m vectors of length p}
#' \item{betam}{beta parameter estimate at final step M}
#' The notation for the return values is described in Ehrlinger and Ishwaran (2012).
#' 
#' @seealso \code{\link{print.l2boost}}, \code{\link{plot.l2boost}}, \code{\link{predict.l2boost}}, 
#' \code{\link{coef.l2boost}}, \code{\link{residuals.l2boost}}, \code{\link{fitted.l2boost}} methods of l2boost 
#' and \code{\link{cv.l2boost}} for K fold cross-validation of the l2boost method. 
#'
#' @examples
#' #--------------------------------------------------------------------------
#' # Example 1: Diabetes data
#' #  
#' # See Efron B., Hastie T., Johnstone I., and Tibshirani R. 
#' # Least angle regression. Ann. Statist., 32:407-499, 2004.
#' data(diabetes, package="l2boost")
#' 
#' l2.object <- l2boost(diabetes$x,diabetes$y, M=1000, nu=.01)
#'
#' # Plot the boosting rho, and regression beta coefficients as a function of
#' # boosting steps m
#' #
#' # Note: The selected coordinate trajectories are colored in red after selection, and 
#' # blue before. Unselected coordinates are colored grey.
#' #
#' par(mfrow=c(2,2))
#' plot(l2.object)
#' plot(l2.object, type="coef")
#' 
#' # increased shrinkage and number of iterations.
#' l2.shrink <- l2boost(diabetes$x,diabetes$y,M=5000, nu=1.e-3) 
#' plot(l2.shrink)
#' plot(l2.shrink, type="coef")
#'
#' \dontrun{
#' #--------------------------------------------------------------------------
#' # Example 2: elasticBoost simulation
#' # Compare l2boost and elastic net boosting
#' # 
#' # See Zou H. and Hastie T. Regularization and variable selection via the 
#' # elastic net. J. Royal Statist. Soc. B, 67(2):301-320, 2005
#' set.seed(1025)
#' 
#' # The default simulation uses 40 covariates with signal concentrated on 
#' # 3 groups of 5 correlated covariates (for 15 signal covariates)
#' dta <- elasticNetSim(n=100)
#' 
#' # l2boost the simulated data with groups of correlated coordinates
#' l2.object <- l2boost(dta$x,dta$y,M=10000, nu=1.e-3, lambda=NULL)
#' 
#' par(mfrow=c(2,2))
#' # plot the l2boost trajectories over all M
#' plot(l2.object, main="l2Boost nu=1.e-3")
#' # Then zoom into the first m=500 steps
#' plot(l2.object, xlim=c(0,500), ylim=c(.25,.5), main="l2Boost nu=1.e-3")
#' 
#' # elasticNet same data with L1 parameter lambda=0.1
#' en.object <- l2boost(dta$x,dta$y,M=10000, nu=1.e-3, lambda=.1) 
#' 
#' # plot the elasticNet trajectories over all M
#' #
#' # Note 2: The elasticBoost selects all coordinates close to the selection boundary,
#' # where l2boost leaves some unselected (in grey)
#' plot(en.object, main="elasticBoost nu=1.e-3, lambda=.1")
#' # Then zoom into the first m=500 steps
#' plot(en.object, xlim=c(0,500), ylim=c(.25,.5),
#'   main="elasticBoost nu=1.e-3, lambda=.1")
#' }
#' 
#' @rdname l2boost
#' @name l2boost
#' @method l2boost default
#' @S3method l2boost default
l2boost.default <- function(x, y,
                            M = NULL, nu = 1e-4, lambda = NULL, trace = FALSE, 
                            type = c("discrete", "hybrid", "friedman","lars"),
                            qr.tolerance = 1e-30,
                            eps.tolerance = .Machine$double.eps,
                            ...) {
  call<-match.call()
  # preliminary checks, set dimensions
  if (nu <= 0 || nu > 1) stop("nu set incorrectly:", nu, "\n")
  x.org <- x <- as.matrix(x)
  y.org <- y <- c(y)
  n.org <- n <- length(y)
  p <- ncol(x)
  if (n < 2)  stop("insufficient data, n =", n, "\n")
  
  # process the data
  ybar <- mean(y.org)
  y <- y - ybar
  cNames <- colnames(x)
  rownames(x) <- colnames(x) <- NULL
  x <- scale(x)/sqrt(n - 1)
  x.attr <- attributes(x)
  
  # quick design matrix sanity check.
  if(all(is.na(x)))
    stop("Completely singular design matrix, probably due to all columns being coincident with each other.")
  
  # remove x-columns with NA standardized values
  x.na <- sapply(1:p, function(j) {any(is.na(x[, j]))})
  if (any(x.na)) {
    x <- x[, !x.na]
    x.org <- x.org[, !x.na]
    p <- ncol(x)
    x <- scale(x)/sqrt(n - 1)
    x.attr <- attributes(x)
  }
  
  # elastic net modification to x, y
  enet <- (!is.null(lambda) && lambda > 0) 
  if (enet) {
    if (lambda < 0) stop("lambda must be positive")
    x <- rbind(x, diag(sqrt(lambda), p))/sqrt(1 + lambda)
    y <- c(y, rep(0, p))
    n <- length(y)
  }
  
  # decide which algorithm is to be applied
  type <- match.arg(type)
  TYPE <- switch(type,
                 friedman = "FRIEDMAN",
                 discrete = "DISCRETE",
                 hybrid   = "HYBRID",
                 lars     = "LARS")
  
  if (trace) cat(paste("implementing l2boost", TYPE, "..."), "\n")
  
  # special scenario for LARS
  if (TYPE == "LARS") {
    if (is.null(M)) M <- min(n, p - 1) else M <- min(M, p - 1)
  }
  else {
    if (is.null(M)) M <- min(n, p)
  }
  
  # initialize predictor
  # initialize beta 
  Fm <- rep(0, n)
  betam <- rep(0, p)
  Fm.path <- betam.path <- vector("list", length = (M + 1))
  Fm.path[[1]] <- Fm
  betam.path[[1]] <- betam
  
  # initialize (l,L,S) triplets
  # initialize rho entities
  # initialize x-correlation (corr) entities
  l.crit <- L.crit <- S.crit <- rep(0, M + 1)
  rhom.path <- vector(length = (M + 1), "list")
  corr.x <- vector(length = p, "list")
  VR <- rep(0, M)
  mjk <- vector(length = M, "list")
  stepSize <- vector(length = M, "list")
  
  # calculate all x_j^Ty: initializes rho
  # initialize rho-related quantities
  rho.m <- c(t(x) %*% y)
  l1 <- l.crit[1] <- resample(which.max.ind(abs(rho.m)))
  corr.x[[l1]] <- extract.corr(x, l1, enet, n.org)
  rhom.path[[1]] <- rho.m
  
  # cycling busting could terminate prior to the final r=M step
  early.terminate.flag <- FALSE
  
  #-----------------------------------------------------
  # Main entry point
  
  for (r in 1:M) {
    
    # assign lr
    lr <- l.crit[r]
    
    # extract the R_j,l correlation: only need new values
    if (r > 1 && (sum(lr == l.crit[1:(r - 1)]) == 0)) {
      corr.x[[lr]] <- extract.corr(x, lr, enet, n.org)
    }
    
    # What we do now depends upon TYPE 
    
    if (TYPE == "FRIEDMAN") {
      
      ##---------------FRIEDMAN-------------------
      ## original, bare-bones l2boost
      
      # get the next critical value
      # update the gradient
      # break ties randomly
      rho.update <- nu * rho.m[lr]
      rho.m.PlusOne <- rho.m - rho.update * corr.x[[lr]]
      lr.PlusOne <- resample(which.max.ind(abs(rho.m.PlusOne)))
      
      # trace
      if (trace) cat(r, lr, lr.PlusOne, nu, "\n")
      
      # updates
      L.crit[r + 1] <- 1
      l.crit[r + 1] <- lr.PlusOne 
      Fm <- Fm + rho.update * x[, lr]
      betam[lr] <- betam[lr] + rho.update
      rho.m <- rho.m.PlusOne
      
      #save mjk (methodology related)
      mjk[[r]] <- mstep.long(rho.m, corr.x[[lr]], lr, nu)
      
    }
    
    else if (TYPE == "DISCRETE") {
      
      ##---------------DISCRETE SOLUTION PATH-------------------
      ## unlike HYBRID:
      #  if nu is small you get one big step, followed by steps of size 1
      ## if nu is large you get steps of size 1
      
      # get the discrete path solution (uses random tie-breaking)
      path.solve <- get.discrete.solution(rho.m, corr.x, lr, nu)
      
      # (l,L,S) update; Vr
      l.crit[r + 1] <- path.solve$lr.PlusOne 
      L.crit[r + 1] <- path.solve$Lr 
      S.crit[r + 1] <- S.crit[r] + path.solve$Lr
      Vr <- (1 - (1 - nu)^path.solve$Lr)
      
      # trace
      if (trace) cat(r, lr, path.solve$lr.PlusOne, Vr, "\n")
      
      # updates
      rho.update <- Vr * rho.m[lr]
      Fm <- Fm + rho.update * x[, lr]
      betam[lr] <- betam[lr] + rho.update
      rho.m <- rho.m - rho.update * corr.x[[lr]]
      
      #save mjk (methodology related)
      mjk[[r]] <- path.solve$M.step
      
    }
    
    else if (TYPE == "HYBRID") {
      
      ##---------------HYBRID-------------------
      ## takes a big first step size even if nu is small
      ## ... works best if nu is moderate and not too small
      
      # get the dynamic path solution
      # breaks ties randomly using the gradient
      path.solve <- get.hybrid.solution(rho.m, corr.x, lr)
      
      # deal with ties
      
      #case 1, Lr>1 but dynamic ties (EXTREMELY RARE)
      if ((path.solve$nu.r.max > nu) & (length(path.solve$lr.PlusOne) > 1)) {
        cat("numerical issue with DYNAMIC ties: case1\n")                
        M.step <- floor(1 + log(1 - path.solve$nu.r)/log(1 - nu))
        lr.PlusOne <-  which.min.ind(M.step, lr)
        Lr <- unique(M.step[lr.PlusOne])
        lr.PlusOne <- break.ties(rho.m, corr.x[[lr]], lr, lr.PlusOne, Lr, nu)
      }
      #case 2, Lr=1, dynamic ties (RARE)
      else if ((path.solve$nu.r.max <= nu) & (length(path.solve$lr.PlusOne) > 1)) {
        cat("numerical issue with DYNAMIC ties: case2\n")
        lr.PlusOne <- break.ties(rho.m, corr.x[[lr]], lr, path.solve$lr.PlusOne, 1, nu)
      }
      else {
        #case 3, everything else (99.99999....% of the time)
        #(do nothing)
        lr.PlusOne <- path.solve$lr.PlusOne
      }
      
      # floor-bound on VR
      VR[r] <- max(path.solve$nu.r[lr.PlusOne], nu, na.rm = TRUE)
      
      # trace
      if (trace) cat(r, lr, lr.PlusOne, path.solve$nu.r[lr.PlusOne], VR[r], "\n")
      
      # updates
      l.crit[r + 1] <- lr.PlusOne 
      rho.update <- VR[r] * rho.m[lr]
      Fm <- Fm + rho.update * x[, lr]
      betam[lr] <- betam[lr] + rho.update
      rho.m <- rho.m - rho.update * corr.x[[lr]]
      
    }
    
    else if (TYPE == "LARS") {
      
      ##---------------LARS-l2boost-LIMIT-------------------
      # get the l2boost-lars-limit
      # first need to define the active set 
      # break step size ties randomly
      active.set <- l.crit[1:r]
      first.coord <- active.set[1]
      path.solve <- get.lars.solution(rho.m, corr.x, active.set, qr.tolerance, eps.tolerance)
      active.set.sign <- path.solve$active.set.sign
      VR[r] <- path.solve$nu.limit
      
      # trace
      if (trace) cat(r, lr, path.solve$lr.PlusOne, VR[r], "\n")
      
      # update step size
      stepSize[[r]] <- abs(path.solve$gmma)/sum(abs(path.solve$gmma), na.rm = TRUE)
      
      # updates
      l.crit[r + 1] <- path.solve$lr.PlusOne 
      rho.update <- VR[r] * rho.m[first.coord]
      Fm <- Fm + rho.update * rowSums(t(t(as.matrix(x[, active.set])) * active.set.sign * path.solve$gmma))
      betam[active.set] <- betam[active.set] + rho.update * active.set.sign * path.solve$gmma
      rho.m <- rho.m - rho.update * rowSums(sapply(1:length(active.set), function(j) {
        corr.x[[active.set[j]]] * path.solve$gmma[j] * active.set.sign[j]}))
      
    }
    
    
    # -----------common path updates---------------------
    # update Fm, betam and rhom paths (add ybar to Fm's path)
    Fm.path[[r + 1]] <- ybar + Fm 
    betam.path[[r + 1]] <- betam
    rhom.path[[r + 1]] <- rho.m
    
  }
  
  #---------check to see if early termination occurred during cycling busting----------
  if (early.terminate.flag) M <- r - 1
  
  #---------elastic net modification----------
  if (enet) {
    n.org <- length(y.org)
    Fm <- Fm[1:n.org]
    Fm.path <- lapply(1:(M + 1), function(r) {Fm.path[[r]][1:n.org]})
    betam <- betam * sqrt(1 + lambda)
    betam.path <- lapply(1:(M + 1), function(r) {betam.path[[r]] * sqrt(1 + lambda)})
  }    
  
  # ---------return the goodies------------------
  # trim M+1 indices
  # adjust Fm by adding ybar
  object <- list(
    call = call, type=type,
    l.crit = l.crit[1:M],
    L.crit = (if (TYPE == "FRIEDMAN" | TYPE == "DISCRETE") L.crit[-1] else VR[1:M]),
    S.crit = (if (TYPE == "DISCRETE") S.crit[-1] else NULL),
    mjk = (if (TYPE == "FRIEDMAN" | TYPE == "DISCRETE") mjk  else NULL),
    stepSize = (if (TYPE == "LARS") stepSize else NULL),        
    rhom.path = rhom.path[1:(M+1)],
    Fm = (ybar + Fm), Fm.path = Fm.path[1:(M+1)],
    betam = betam, betam.path = betam.path[1:(M+1)],
    x = x.org, x.attr = x.attr, x.na = x.na, names = cNames,
    y = y.org, ybar = ybar, nu=nu, lambda=lambda)
  class(object) <- c("l2boost", TYPE)
  invisible(object)
  
}

#' 
#' @usage \method{l2boost}{formula}(formula, data, ...)
#' @name l2boost
#' @rdname l2boost
#' 
#' @aliases l2boost.formula l2boost.default
#' 
#' @method l2boost formula
#' @S3method l2boost formula
l2boost.formula <- function(formula, data, ...){  
    mf <- model.frame(formula=formula, data=data)
    x<- model.matrix(attr(mf, "terms"), data=mf)
    y<-model.response(mf)
    
    est<- l2boost.default(x,y, ...)
    est$call<-match.call()
    est$formula <-formula
    invisible(est)
  }
