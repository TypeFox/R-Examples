##' Compute (on a regular grid) the empirical attainment function from conditional simulations 
##' of Gaussian processes corresponding to two objectives. This is used to estimate the Vorob'ev
##' expectation of the attained set and the Vorob'ev deviation.
##'
##' @title Conditional Pareto Front simulations
##' 
##' @param fun1sims numeric matrix containing the conditional simulations of the first output (one sample in each row),
##' @param fun2sims numeric matrix containing the conditional simulations of the second output (one sample in each row),
##' @param response a matrix containing the value of the two objective functions, one
##' output per row,
##' @param paretoFront optional matrix corresponding to the Pareto front of the observations. It is 
##' estimated from \code{response} if not provided,
##' @param f1lim optional vector (see details),
##' @param f2lim optional vector (see details),
##' @param refPoint optional vector (see details),
##' @param n.grid integer determining the grid resolution,
##' @param compute.VorobExp optional boolean indicating whether the Vorob'ev Expectation
##' should be computed. Default is \code{TRUE},
##' @param compute.VorobDev optional boolean indicating whether the Vorob'ev deviation
##' should be computed. Default is \code{TRUE}.
##' @return A list which is given the S3 class "\code{CPF}".
##' \itemize{
##'  \item{\code{x, y}}{: locations of grid lines at which the values of the attainment
##'   are computed,}
##'  \item{\code{values}}{: numeric matrix containing the values of the attainment on the grid,}
##'  \item{\code{PF}}{: matrix corresponding to the Pareto front of the observations,}
##'  \item{\code{responses}}{: matrix containing the value of the two objective functions, one
##' objective per column,}
##' \item{\code{fun1sims, fun2sims}}{: conditional simulations of the first/second output,}
##' \item{\code{VE}}{: Vorob'ev expectation, computed if \code{compute.VorobExp = TRUE} (default),}
##' \item{\code{beta_star}}{: Vorob'ev threshold, computed if \code{compute.VorobExp = TRUE} (default),}
##' \item{\code{VD}}{: Vorov'ev deviation, computed if \code{compute.VorobDev = TRUE} (default),}
##' }
##' @details Works with two objectives. The user can provide locations of grid lines for
##'  computation of the attainement function with vectors \code{f1lim} and \code{f2lim}, in the form of regularly spaced points. 
##'  It is possible to provide only \code{refPoint} as a reference for hypervolume computations.
##'  When missing, values are determined from the axis-wise extrema of the simulations.
##' @seealso Methods \code{coef}, \code{summary} and \code{plot} can be used to get the coefficients from a \code{CPF} object, 
##' to obtain a summary or to display the attainment function (with the Vorob'ev expectation if \code{compute.VorobExp} is \code{TRUE}). 
##' @importFrom emoa nondominated_points is_dominated dominated_hypervolume
##' @import DiceKriging
##' @importFrom rgenoud genoud
##' @importFrom grDevices gray
##' @importFrom graphics filled.contour axis points contour lines polygon par
##' @references
##' M. Binois, D. Ginsbourger and O. Roustant (2015), Quantifying Uncertainty on Pareto Fronts with Gaussian process conditional simulations, 
##' \emph{European Journal of Operational Research}, 243(2), 386-394. \cr \cr
##' C. Chevalier (2013), \emph{Fast uncertainty reduction strategies relying on Gaussian process models}, University of Bern, PhD thesis. \cr \cr
##' I. Molchanov (2005), \emph{Theory of random sets}, Springer.
##' @export
##' @examples
##' library(DiceDesign)
##' set.seed(42)
##' 
##' nvar <- 2
##' 
##' fname <- "P1" # Test function
##' 
##' # Initial design
##' nappr <- 10
##' design.grid <- maximinESE_LHS(lhsDesign(nappr, nvar, seed = 42)$design)$design
##' response.grid <- t(apply(design.grid, 1, fname))
##' 
##' # kriging models: matern5_2 covariance structure, linear trend, no nugget effect
##' mf1 <- km(~., design = design.grid, response = response.grid[,1])
##' mf2 <- km(~., design = design.grid, response = response.grid[,2])
##' 
##' # Conditional simulations generation with random sampling points 
##' nsim <- 40
##' npointssim <- 150 # increase for better results
##' Simu_f1 <- matrix(0, nrow = nsim, ncol = npointssim)
##' Simu_f2 <- matrix(0, nrow = nsim, ncol = npointssim)
##' design.sim <- array(0, dim = c(npointssim, nvar, nsim))
##' 
##' for(i in 1:nsim){
##'   design.sim[,,i] <- matrix(runif(nvar*npointssim), nrow = npointssim, ncol = nvar)
##'   Simu_f1[i,] <- simulate(mf1, nsim = 1, newdata = design.sim[,,i], cond = TRUE,
##'                           checkNames = FALSE, nugget.sim = 10^-8)
##'   Simu_f2[i,] <- simulate(mf2, nsim = 1, newdata = design.sim[,,i], cond = TRUE, 
##'                           checkNames = FALSE, nugget.sim = 10^-8)
##' }
##' 
##' # Attainment and Voreb'ev expectation + deviation estimation
##' CPF1 <- CPF(Simu_f1, Simu_f2, response.grid)
##' 
##' # Details about the Vorob'ev threshold  and Vorob'ev deviation
##' summary(CPF1)
##' 
##' # Graphics
##' plot(CPF1)
##' 

CPF <-
  function(fun1sims, fun2sims, response, paretoFront = NULL, f1lim = NULL, f2lim = NULL, refPoint = NULL,
           n.grid = 100, compute.VorobExp = TRUE, compute.VorobDev = TRUE){
    
    ## @param compute.VorobExpect an optional boolean indicating whether the Vorob'ev expectation
    ## should be computed after the attainment , compute.VorobExpect = TRUE
    
    if(is.null(paretoFront)){
      paretoFront <- t(nondominated_points(t(response)))
    }
    
    # Determination of the grid for computation of the attainment (if X is NULL)
    ##1) refPoint is provided : take the lower bound from the simulations
    ##2) Nothing is provided : take the bounds from the simulations
    
    if(is.null(refPoint)){
      if(is.null(f1lim)){
        f1lim <- seq(min(fun1sims), max(fun1sims), length.out = n.grid)
      }
      if(is.null(f2lim)){
        f2lim <- seq(min(fun2sims), max(fun2sims), length.out = n.grid)
      }
    }else{
      if(is.null(f1lim)){
        f1lim <- seq(min(fun1sims), refPoint[1], length.out = n.grid)
      }
      if(is.null(f2lim)){
        f2lim <- seq(min(fun1sims), refPoint[2], length.out = n.grid)
      }
    }
    
    X <- as.matrix(expand.grid(f1lim, f2lim))
    
    nsim <- dim(fun1sims)[1]
    
    Pdom <- rep(0,dim(X)[1])
    
    ## Some points on the grid can be computed quickly
    ActivePoints <- rep(0,dim(X)[1])
    
    ## for the points which dominates the non-dominated points of each simulation, the probability is 0 
    ActivePoints[which(X[,1] <= min(fun1sims) | X[,2] <= min(fun2sims))] = 1
    
    ## for the points in X which are dominated by the Pareto Front of the observations, the probability is 1
    for(i in 1:dim(paretoFront)[2]){
      tmp <- which(X[,1] >= paretoFront[i,1] & X[,2] >= paretoFront[i,2])
      Pdom[tmp] = nsim
      ActivePoints[tmp] = 1
    }
    
    
    for(i in 1:nsim){
      tmp <- nondominated_points(cbind(rbind(fun1sims[i,],fun2sims[i,]), t(paretoFront)))
      for(j in 1:dim(X)[1]){
        if(ActivePoints[j] == 0 && is_dominated(cbind(X[j,], tmp))[1]){
          Pdom[j] = Pdom[j]+1
        }
      }  
    }
    
    res <- list(x = f1lim, y = f2lim, values = matrix(Pdom/nsim, length(f1lim), length(f2lim)),
                PF = paretoFront, responses = response, fun1sims = fun1sims, 
                fun2sims = fun2sims, beta_star = NA, VE = NULL, VD = NA) 
    
    if(compute.VorobExp){
      res$beta_star <- VorobThreshold(res)
      res$VE <- VorobExpect(res)
    }
    
    if(compute.VorobDev){
      if(compute.VorobExp == FALSE){
        stop("Vorob'ev expectation required, set compute.VorobExp to TRUE")
      }
      res$VD <- VorobDev(res)
    }
    
    class(res) <- "CPF"
    res
  }


##' @method summary CPF
##' @export
summary.CPF <- function(object,...){
  ans <- object
  class(ans) <- "summary.CPF"
  ans
}



##' @method plot CPF
##' @export
plot.CPF <- function(x, ...){
  # #' @title Plot of the attainment function and Vorob'ev expectation.
  # #' @details See \code{\link[DiceMOO]{CPF}} for an example.
  
  if(is.null(x$VE)){
    filled.contour(x$x, x$y, x$values,
                   col=graypalette(7), levels = c(0,0.01,0.2,0.4,0.6,0.8,0.99,1),
                   main = "Empirical attainment function",
                   xlab = expression(f[1]), ylab = expression(f[2]),
                   plot.axes = { axis(1); axis(2);
                                 #points(x$PF[1, ], x$PF[2, ], pch = 17, col = "blue");
                                 points(x$response[,1], x$responses[,2], pch = 17, col="blue");
                                 contour(x$x, x$y, x$values, add=T,
                                         levels = c(0,0.2,0.4,0.6,0.8,1),labcex=1);
                                 plotParetoEmp(x$PF,col="blue", lwd = 2);
                   }
    )
  }else{
    filled.contour(x$x, x$y, x$values,
                   col=graypalette(7), levels = c(0,0.01,0.2,0.4,0.6,0.8,0.99,1),
                   main = substitute(paste("Empirical attainment function, ",beta,"* = ", a),
                                     list(a = formatC(x$beta_star, digits = 2))),
                   xlab = expression(f[1]), ylab = expression(f[2]),
                   plot.axes = { axis(1); axis(2);
                                 #points(x$PF[1, ], x$PF[2, ], pch = 17, col = "blue");
                                 points(x$response[,1], x$responses[,2], pch = 17, col="blue");
                                 contour(x$x, x$y, x$values, add=T,
                                         levels = c(0,0.2,0.4,0.6,0.8,1),labcex=1);
                                 plotParetoEmp(x$PF,col="blue", lwd = 2);
                                 plotParetoEmp(x$VE,col="cyan",lwd=2)
                   }
    )
  }
}

#' @export
print.summary.CPF <- function(x, ...) {
  ## Number of simulations
  cat("\n Number of simulations used:", dim(x$fun1sims)[1], "\n")
  
  ## Number of simulations points
  cat("\n Number of simulations points:", dim(x$fun1sims)[2], "\n")
  
  ## Ref(Nadir) Point
  cat("\n Ref Point:", max(x$x), max(x$y), "\n")
  
  ## Ideal Point
  cat("\n Ideal Point:", min(x$x), min(x$y), "\n")
  
  ## beta star
  if(is.na(x$beta_star)){
    cat("\n Vorob'ev threshold not computed")
  }else{
    cat("\n Vorob'ev threshold:", x$beta_star, "\n")
  }
  
  ## Vorob'ev deviation
  if(is.na(x$VD)){
    cat("\n Vorob'ev deviation not computed")
  }else{
    cat("\n Vorob'ev deviation:", x$VD, "\n")
  }
}

##' @method coef CPF
##' @export
coef.CPF <- function(object, type = c("grid", "attainment", "threshold",
                                      "expectation", "refPoint", "deviation"), ...){
  # #' Get coefficients values from a CPF object.
  # #' @title Get coefficients values
  # #' @param object an object with class \code{\link[DiceMOO]{CPF}} from which the coefficiant will be extracted.
  # #' @param type Character string or vector specifying which type(s) of coefficients in the structure will be
  # #'  extracted. Should be one of \code{grid}, \code{attainment}, \code{threshold}, \code{expectation}, \code{refPoint}, \code{deviation}.
  # #' @param ... other arguments (undocumented at this stage). 
  type <- match.arg(type)
  switch(type,
         grid = as.matrix(expand.grid(object$x, object$y)),
         attainment = object$values,
         threshold = object$beta_star,
         expectation = object$VE,
         refPoint = c(max(object$x), max(object$y)),
         deviation = object$VD
  )
}


graypalette <- function(n){
  sapply(seq(0,1,,n)^0.5,gray)
}

# ' Compute the Vorob'ev threshold (on a grid)
# ' @title Compute the Vorob'ev threshold
# ' @param Attainment List obtained using the empAttainFun
# ' @param precision precision on Beta*
# ' @value beta_star Vorob'ev threshold
# ' @details Based on a discrete approximation, the Vorob'ev threshold is obtained from the condition
# ' mu(Q_beta) <= E(mu(X)) <= mu(Q_beta*) forall beta > beta* . 
# ' 
VorobThreshold <- function(A, precision = 0.001){
  
  Q_beta_star <- matrix(0, length(A$x), length(A$y))
  Q_beta_star[which(A$values > 0.5)] = 1
  
  a <- 0
  b <- 1
  
  # Dichotomy
  beta_star <- 0.5
  tmp <- sum(Q_beta_star)
  
  # The integral of the coverage function gives the expected volume
  Emu <- sum(A$values)
  
  # The dichotomy stops when a plateau is reached or the difference between to 
  
  while(b - a > precision){
    if(tmp == Emu) 
      break
    
    if(tmp < Emu){
      b <- (a+b)/2
    }else{
      a <- (a+b)/2
    }
    beta_star = (a+b)/2
    Q_beta_star <- matrix(0, length(A$x), length(A$y))
    Q_beta_star[which(A$values > beta_star)] = 1
    tmp = sum(Q_beta_star)
  }
  return(beta_star)
}


# ' Compute the Vorob'ev expectation
# ' @title Compute the Vorob'ev expectation
# ' @param List obtained using the empAttainFun
# ' @param beta_star Vorob'ev threshold
# ' @value VorobevExpect Pareto frontier of the Vorob'ev expectation

VorobExpect <- function(A){
  
  VE <- rep(0,length(A$y))
  
  for(i in length(A$y):1){
    ind <- which(A$values[,i] < A$beta_star)
    if(length(ind) == 0){ # case when the box is too small
      VE[i] <- NA
    }else{
      if(max(ind) == length(A$x) && A$values[length(A$x), i + 1] <= A$beta_star){
        VE[i] <- NA
      }else{
        VE[i] <- A$x[max(ind)]
      }
    }
  }
  
  VE <- cbind(VE, A$y)
  
  VE <- t(nondominated_points(t(VE[which(!is.na(VE[,1])),])))
  
  return(VE)
}


# # ' Compute the Vorob'ev deviation of Random Non-dominated Points sets
# # ' @title Vorob'ev deviation of RNP sets
# # ' @details
# # ' In this case there is no need for a grid for computation of volumes : only the hypervolume is necessary
# # ' @examples
# # ' #------------------------------------------------------------
# # ' # Simple example
# # ' #------------------------------------------------------------
# # ' VE <- rbind(c(1,3,5), c(3,2,1))
# # ' f1 <- rbind(c(2,4,6), c(0,1,1), c(0,1,4))
# # ' f2 <- rbind(c(3,2,1), c(3,2,1), c(3,2,4))
# # ' 
# # ' VorobDev(f1, f2, VE, refPoint = c(7,4))
# # ' 
# # ' # Graphical verification
# # ' plot(t(VE), type = "p", xlim = c(-1,8), ylim = c(-1,5),
# # '  xlab = expression(f[1]), ylab = expression(f[2]),
# # '  pch = 20, col ="cyan", asp = 1)
# # ' plotParetoEmp(VE, col = "cyan", lwd = 2)
# # ' 
# # ' points(7,4, pch = 3, lwd = 4)
# # ' 
# # ' points(f1[1,], f2[1,], pch = 20 , col = "red")
# # ' plotParetoEmp(nondominated_points(t(cbind(f1[1,],f2[1,]))), col = "red")
# # ' points(f1[2,], f2[2,], pch = 20 , col = "green")
# # ' plotParetoEmp(nondominated_points(t(cbind(f1[2,],f2[2,]))), col = "green")
# # ' points(f1[3,], f2[3,], pch = 20 , col = "orange")
# # ' plotParetoEmp(nondominated_points(t(cbind(f1[3,],f2[3,]))), col = "orange")
# # ' @export
VorobDev <- function(CPF, refPoint = NULL){
  
  nsim = dim(CPF$fun1sims)[1]
  
  # Determination of the grid for computation of the attainment (if X is NULL)
  ##1) refPoint is provided : take the lower bound from the simulations
  ##2) Nothing is provided : take the bounds from the simulations
  
  if(is.null(refPoint)){
    refPoint <- c(max(CPF$fun1sims), max(CPF$fun2sims))
  }
  
  mi1 <- min(CPF$fun1sims)
  ma1 <- refPoint[1]
  
  mi2 <- min(CPF$fun2sims)
  ma2 <- refPoint[2]
  
  #lowPoint <- c(min(fun1sims), min(fun2sims))
  
  VD <- 0 #Vorob'ev deviation
  
  # symmetric difference : 2*H(AUB) - H(A) - H(B)
  for(i in 1:nsim){
    H1 <- dominated_hypervolume(rbind(CPF$fun1sims[i,],CPF$fun2sims[i,]), ref = c(ma1, ma2))
    H2 <- dominated_hypervolume(t(CPF$VE), ref = c(ma1, ma2))
    H12 <- dominated_hypervolume(cbind(rbind(CPF$fun1sims[i,],CPF$fun2sims[i,]),
                                       t(CPF$VE)), ref = c(ma1, ma2))
    VD = VD + 2*H12 -H1 -H2
    
  }
  
  return(VD/nsim)
}
