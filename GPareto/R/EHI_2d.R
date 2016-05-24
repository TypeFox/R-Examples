## ' Computes the analytical Expected Hypervolume Improvement with respect to the
## ' current Pareto front. To avoid numerical instabilities, the new point is evaluated only if it is not too close to an existing observation.
## ' 
## ' @title Analytical expression of the Expected Hypervolume Improvement with 2 objectives
## ' 
## ' @param x a vector representing the input for which one wishes to calculate EHI,
## ' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
## ' @param critcontrol list with element \code{refPoint} (reference point) for hypervolume computations.
## ' Options for the \code{\link[GPareto]{checkPredict}} function: \code{threshold} and \code{distance} are available to avoid numerical issues occuring when adding points too close to the existing ones.
## ' @param type "SK" or "UK" (by default), depending whether uncertainty related to trend estimation 
## '        has to be taken into account. 
## ' @param paretoFront matrix corresponding to the Pareto Front (one output per column).
## ' @return The expected hypervolume improvement at \bold{x}.
## ' @seealso \code{\link[DiceOptim]{EI}} from package DiceOptim from which \code{EHI} is an extension 
## ' @details Both models are supposed to share the same design. \cr
## '          It is adapted from the Matlab source code by Michael Emmerich and Andre Deutz, LIACS, 
## '          Leiden University, 2010 available here :
## '          \url{http://natcomp.liacs.nl/code/HV_based_expected_improvement.zip}.
## ' @export
## ' @useDynLib GPareto
## ' @importFrom Rcpp evalCpp
## ' @references Emmerich, M., & Klinkenberg, J. W. (2008). The computation of the expected improvement 
## ' in dominated hypervolume of Pareto front approximations. Leiden Institute for Advanced Computer Science, 
## ' Tech. Rep, 1.
## ' @examples
## ' #---------------------------------------------------------------------------
## ' # EHI surface associated with the "P1" problem at a 15 points design
## ' #---------------------------------------------------------------------------
## ' \donttest{
## ' set.seed(25468)
## ' library(DiceDesign)
## ' 
## ' n_var <- 2 
## ' f_name <- "P1" 
## ' n.grid <- 21
## ' test.grid <- expand.grid(seq(0, 1,, n.grid), seq(0, 1,, n.grid))
## ' n_appr <- 15 
## ' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
## ' response.grid <- t(apply(design.grid, 1, f_name))
## ' Front_Pareto <- t(nondominated_points(t(response.grid)))
## ' mf1 <- km(~., design = design.grid, response = response.grid[,1])
## ' mf2 <- km(~., design = design.grid, response = response.grid[,2])
## ' 
## ' EHI_grid <- apply(test.grid, 1, EHI_2d, model = list(mf1, mf2),
## '                   critcontrol = list(refPoint = c(300,0)))
## ' 
## ' filled.contour(seq(0, 1,, n.grid), seq(0, 1,, n.grid), matrix(EHI_grid, n.grid),
## '                main = "Expected Hypervolume Improvement", xlab = expression(x[1]),
## '                ylab = expression(x[2]), color = terrain.colors, nlevels = 50,
## '                plot.axes = {axis(1); axis(2);
## '                             points(design.grid[,1],design.grid[,2],pch=21,bg="white")
## '                             }
## '               )
## ' }

EHI_2d <- function(x, model, critcontrol=NULL, type = "UK", paretoFront = NULL){
  
  n.obj <- length(model)
  d <- model[[1]]@d
  x.new <- matrix(x, 1, d)
  
  if(is.null(paretoFront) || is.null(critcontrol$refPoint)){
    observations <- matrix(0, model[[1]]@n, n.obj)
    for (i in 1:n.obj) observations[,i] <- model[[i]]@y
    if(is.null(paretoFront))
      paretoFront <- t(nondominated_points(t(observations)))
  }
  
  if (is.unsorted(paretoFront[,1])){
    paretoFront <- paretoFront[order(paretoFront[,1]),]
  }
  
  refPoint <- critcontrol$refPoint
  if (is.null(refPoint)){
    refPoint <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj) ### May be changed? 
    cat("No refPoint provided, ", signif(refPoint, 3), "used \n")
  } 
  
  if (n.obj!=2){
    print("Analytical hypervolume EI only works with 2 objectives")
    return(NULL)
  } else {
    
    mu    <- rep(NaN, n.obj)
    sigma <- rep(NaN, n.obj)
    for (i in 1:n.obj){    
      pred     <- predict(object=model[[i]], newdata=x.new, type=type, checkNames = FALSE, light.return = TRUE)
      mu[i]    <- pred$mean
      sigma[i] <- pred$sd
    }
    
    ## A new x too close to the known observations could result in numerical problems
  
    if(checkPredict(x, model, type = type, distance = critcontrol$distance, threshold = critcontrol$threshold)){
      return(-1)
    }else{
      return(EHI_2d_wrap_Rcpp(paretoFront, refPoint, mu, sigma))
    }
  }
}
