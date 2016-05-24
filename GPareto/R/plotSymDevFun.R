##' Display the Symmetric Deviation Function for an object of class CPF.
##' @title Display the Symmetric Deviation Function
##' @param CPF CPF object, see \code{\link[GPareto]{CPF}},
##' @param n.grid number of divisions of the grid in each dimension.
##' @details Display observations in red and the corresponding Pareto front by a step-line. 
##' The blue line is the estimation of the location of the Pareto front of the kriging models, named Vorob'ev expectation. 
##' In grayscale is the intensity of the deviation (symmetrical difference) from the Vorob'ev expectation for couples of conditional simulations.
##' @export
##' @references
##' M. Binois, D. Ginsbourger and O. Roustant (2015), Quantifying Uncertainty on Pareto Fronts with Gaussian process conditional simulations, 
##' \emph{European Journal of Operational Research}, 243(2), 386-394. \cr \cr
##' @examples
##' library(DiceDesign)
##' set.seed(42)
##' 
##' nvar <- 2
##' 
##' # Test function
##' fname = "P1"
##' 
##' # Initial design
##' nappr <- 10
##' design.grid <- maximinESE_LHS(lhsDesign(nappr, nvar, seed = 42)$design)$design
##' response.grid <- t(apply(design.grid, 1, fname))
##' 
##' ParetoFront <- t(nondominated_points(t(response.grid)))
##' 
##' # kriging models : matern5_2 covariance structure, linear trend, no nugget effect
##' mf1 <- km(~., design = design.grid, response = response.grid[, 1])
##' mf2 <- km(~., design = design.grid, response = response.grid[, 2])
##' 
##' # Conditional simulations generation with random sampling points 
##' nsim <- 10 # increase for better results
##' npointssim <- 80 # increase for better results
##' Simu_f1 = matrix(0, nrow = nsim, ncol = npointssim)
##' Simu_f2 = matrix(0, nrow = nsim, ncol = npointssim)
##' design.sim = array(0,dim = c(npointssim, nvar, nsim))
##' 
##' for(i in 1:nsim){
##'   design.sim[,, i] <- matrix(runif(nvar*npointssim), npointssim, nvar)
##'   Simu_f1[i,] = simulate(mf1, nsim = 1, newdata = design.sim[,, i], cond = TRUE,
##'                          checkNames = FALSE, nugget.sim = 10^-8)
##'   Simu_f2[i,] = simulate(mf2, nsim = 1, newdata = design.sim[,, i], cond=TRUE, 
##'                          checkNames = FALSE, nugget.sim = 10^-8)
##' }
##' 
##' # Attainment, Voreb'ev expectation and deviation estimation
##' CPF1 <- CPF(Simu_f1, Simu_f2, response.grid, ParetoFront)
##'
##' # Symmetric deviation function
##' plotSymDevFun(CPF1)
##' 
plotSymDevFun <- function(CPF, n.grid = 100){
  
  flim1 <- range(CPF$x)
  flim2 <- range(CPF$y)
  
  SymDevGraph = matrix(0,n.grid, n.grid)
  for(i in 1:dim(CPF$fun1sims)[1]){
    SymDevGraph = SymDevGraph +
      SymDiffRaw(nondominated_points(rbind(CPF$fun1sims[i,], CPF$fun2sims[i,])), 
                 nondominated_points(t(CPF$VE)),c(flim1[1],flim1[2]),
                 c(flim2[1],flim2[2]),n.grid)
  }
  filled.contour(seq(flim1[1], flim1[2],,n.grid), seq(flim2[1], flim2[2],,n.grid),
                 matrix(SymDevGraph/dim(CPF$fun1sims)[1],n.grid,n.grid),
                 col = gray(11:0/11), levels=seq(0,1,,11),
                 main ="Symmetric deviation function",
                 xlab=expression(f[1]),ylab=expression(f[2]),
                 
                 plot.axes={axis(1);axis(2);  
                            points(CPF$responses[,1], CPF$responses[,2], pch = 17, 
                                   col = "red")
                            plotParetoEmp(CPF$VE, col="blue", lwd=2)
                            plotParetoEmp(CPF$PF, col="red", lwd=2);
                 }
  )
}

# Value is 1 for points in grid which are in the symmetric difference of set1 and set2
SymDiffRaw <- function(set1,set2, xlim, ylim, resGrid){
  
  grille <- as.matrix(expand.grid(seq(xlim[1],xlim[2],,resGrid),seq(ylim[1],ylim[2],,resGrid)))
  
  matRes <- rep(0,resGrid*resGrid)
  
  for (i in 1:dim(grille)[1]){
    a = is_dominated(cbind(grille[i,], set1))[1]
    b = is_dominated(cbind(grille[i,], set2))[1]
    if(a != b)
      matRes[i] = 1
  }
  
  
  return(matRes)
}

