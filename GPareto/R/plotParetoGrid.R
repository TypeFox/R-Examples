#' Plot the Pareto front and set for 2 variables 2 objectives test problems with evaluations on a grid.
#' 
#' @title Visualisation of Pareto front and set
#' 
#' @param fname name of the function considered,
#' @param xlim,ylim numeric vectors of length 2, giving the \code{x} and \code{y} coordinates ranges, default is \code{[0,1] x [0,1]},
#' @param n.grid number of divisions of the grid in each dimension.
#' @export 
#' @examples
#' #------------------------------------------------------------
#' # Examples with test functions
#' #------------------------------------------------------------
#' 
#' plotParetoGrid("ZDT3", n.grid = 21)
#' 
#' plotParetoGrid("P1", n.grid = 21)
#' 
#' plotParetoGrid("MOP2", xlim = c(0, 1), ylim = c(0, 1), n.grid = 21) 



plotParetoGrid <- function(fname = "ZDT1", xlim = c(0,1), ylim = c(0,1), n.grid = 100){
  Plan <- expand.grid(seq(xlim[1], xlim[2], length.out = n.grid),
                      seq(ylim[1], ylim[2], length.out = n.grid))
  
  response.grid <- apply(Plan, 1, fname)
  
  ParetoFront <- t(nondominated_points(response.grid))
  response.grid <- t(response.grid)
  
  ParetoSet <- Plan[which(response.grid[,1] %in% ParetoFront[,1] & response.grid[,2] %in% ParetoFront[,2]),]
  
  par(mfrow = c(1,2))
  
  plot(response.grid, pch = '.', xlab = expression(f[1]), ylab = expression(f[2]), main = "Pareto Front")
  points(ParetoFront, col = "red", pch = 20) 
  plotParetoEmp(ParetoFront, col = "red", lty = 2) 
  
  plot(Plan,pch = '.', xlab = expression(x[1]), ylab = expression(x[2]), main = "Pareto Set")
  points(ParetoSet, col = "red", pch = 20)
  
  par(mfrow = c(1,1))
}