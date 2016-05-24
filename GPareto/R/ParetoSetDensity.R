##' Estimation of the density of Pareto optimal points in the variable space.
##' @title Estimation of Pareto set density
##' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
##' @param lower vector of lower bounds for the variables,
##' @param upper vector of upper bounds for the variables,
##' @param CPS optional matrix containing points from Conditional Pareto Set Simulations 
##' (in the variable space), see details
##' @param nsim optional number of conditional simulations to perform if \code{CPS} is not provided,
##' @param simpoints (optional) If \code{CPS} is \code{NULL}, either a number of simulation points,
##' or a matrix where conditional simulations are to be performed. In the first case, 
##' then simulation points are taken as a maximin LHS design using \code{\link[DiceDesign]{lhsDesign}}.    
##' @param ... further arguments to be passed to \code{\link[ks]{kde}}. 
##' In particular, if the input dimension is greater than three, 
##' a matrix \code{eval.points} can be given (else it is taken as the simulation points). 
##' @return An object of class \code{\link[ks]{kde}} accounting for the 
##' estimated density of Pareto optimal points.
##' @details 
##' This functino estimates the density of Pareto optimal points in the variable space
##' given by the surrogate models. Based on conditional simulations of the objectives at simulation points,
##' Conditional Pareto Set (CPS) simulations are obtained, out of which a density is fitted. 
##'  
##' This function relies on the \code{\link[ks]{ks}} package for the kernel density estimation.
##' @importFrom ks kde
##' @importFrom emoa nds_rank
##' @export
##' @examples
##' \dontrun{ 
##' #---------------------------------------------------------------------------
##' # Example of estimation of the density of Pareto optimal points
##' #---------------------------------------------------------------------------
##' set.seed(42)
##' n_var <- 2 
##' fname <- P1
##' lower <- rep(0, n_var)
##' upper <- rep(1, n_var)
##' 
##' res1 <- easyGParetoptim(fn = fname, lower = lower, upper = upper, budget = 15, 
##' control=list(method = "EHI", inneroptim = "pso", maxit = 20))
##' 
##' estDens <- ParetoSetDensity(res1$history$model, lower = lower, upper = upper)
##' 
##' # graphics
##' par(mfrow = c(1,2))
##' plot(estDens, display = "persp", xlab = "X1", ylab = "X2")
##' plot(estDens, display = "filled.contour2", main = "Estimated density of Pareto optimal point")
##' points(res1$history$model[[1]]@@X[,1], res1$history$model[[2]]@@X[,2], col="blue")
##' points(estDens$x[, 1], estDens$x[, 2], pch = 20, col = rgb(0, 0, 0, 0.15))
##' par(mfrow = c(1,1))
##' }
ParetoSetDensity <- function(model, lower, upper, CPS = NULL, nsim = 50, simpoints = 1000, ...){
  # No CPS provided
  if(is.null(CPS)){
    if(is.null(dim(simpoints))){
      simpoints <- lhsDesign(simpoints, model[[1]]@d)$design
    }
    # add known points to simpoints
    simpoints <- rbind(simpoints, model[[1]]@X)
    # exclude duplicates 
    simpoints <- simpoints[!duplicated(simpoints),]
    
    Simus <- array(0, dim = c(length(model), nsim, nrow(simpoints)))
    
    for(i in 1:length(model)){
      Simus[i,,] <- simulate(model[[i]], nsim = nsim, newdata = simpoints, cond = TRUE,
                             checkNames = FALSE, nugget.sim = 10^-8)
    }
    
    # Creation of CPS
    for (i in 1:nsim){
      non_dom_order <- nds_rank(Simus[,i,])
      CPS <- rbind(simpoints[which(non_dom_order == 1),], CPS)
    }
   }
  if(model[[1]]@d > 3 & !exists("eval.points")){
    fhat <- kde(x = CPS, xmin = lower, xmax = upper, eval.points = simpoints, ...)
  }else{
    fhat <- kde(x = CPS, xmin = lower, xmax = upper, ...)
  }
  
  return(fhat)
} 