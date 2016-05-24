##' Multi-objective optimization and quantification of uncertainty on Pareto fronts, using Gaussian process models.
##' @title Package GPareto
##' @author Mickael Binois, Victor Picheny
##' @docType package
##' @name GPareto
##' @references
##' M. Binois, D. Ginsbourger and O. Roustant (2015), Quantifying Uncertainty on Pareto Fronts with Gaussian process conditional simulations, 
##' \emph{European Journal of Operational Research}, 243(2), 386-394. \cr \cr
##' O. Roustant, D. Ginsbourger and Yves Deville (2012), DiceKriging, DiceOptim:
##'  Two R Packages for the Analysis of Computer Experiments by Kriging-Based Metamodeling and Optimization,
##'   \emph{Journal of Statistical Software}, 51(1), 1-55, \url{http://www.jstatsoft.org/v51/i01/}. \cr \cr
##'  M. T. Emmerich, A. H. Deutz, J. W. Klinkenberg (2011), Hypervolume-based expected improvement: Monotonicity properties and exact computation,
##' \emph{Evolutionary Computation (CEC)}, 2147-2154. \cr \cr
##' V. Picheny (2015), Multiobjective optimization using Gaussian process emulators via stepwise uncertainty reduction, 
##' \emph{Statistics and Computing}, 25(6), 1265-1280. \cr \cr
##' T. Wagner, M. Emmerich, A. Deutz, W. Ponweiser (2010), On expected-improvement criteria for model-based multi-objective optimization,   
##' \emph{Parallel Problem Solving from Nature}, 718-727, Springer, Berlin. \cr \cr
##' J. D. Svenson (2011), \emph{Computer Experiments: Multiobjective Optimization and Sensitivity Analysis}, Ohio State University, PhD thesis. \cr \cr
##' C. Chevalier (2013), \emph{Fast uncertainty reduction strategies relying on Gaussian process models}, University of Bern, PhD thesis.  
##' @details 
##' Important functions: \cr
##' \code{\link[GPareto]{GParetoptim}} \cr
##' \code{\link[GPareto]{easyGParetoptim}} \cr
##' \code{\link[GPareto]{crit_optimizer}} \cr
##' \code{\link[GPareto]{plotGPareto}}\cr
##' \code{\link[GPareto]{CPF}}
##' @note 
##' Part of this work has been conducted within the frame of the ReDice Consortium,
##' gathering industrial (CEA, EDF, IFPEN, IRSN, Renault) and academic 
##' (Ecole des Mines de Saint-Etienne, INRIA, and the University of Bern) partners around
##'  advanced methods for Computer Experiments. (http://www.redice-project.org/).\cr
##' 
##' 
##' The authors would like to thank Yves Deville for his precious advices in R programming and packaging, 
##' as well as Olivier Roustant and David Ginsbourger for testing and suggestions of improvements for this package.
##' We would also like to thank Tobias Wagner for providing his Matlab codes for the SMS-EGO strategy.
##' 
##' @seealso \code{\link[DiceKriging]{DiceKriging}}, \code{\link[DiceOptim]{DiceOptim}}
##' @examples
##' \dontrun{
##' #------------------------------------------------------------
##' # Example 1 : Surrogate-based multi-objective Optimization with postprocessing
##' #------------------------------------------------------------
##' set.seed(25468)
##' 
##' d <- 2 
##' fname <- P2
##' 
##' plotParetoGrid(P2) # For comparison
##' 
##' # Optimization
##' budget <- 25 
##' lower <- rep(0, d) 
##' upper <- rep(1, d)
##'      
##' omEGO <- easyGParetoptim(fn = fname, budget = budget, lower = lower, upper = upper)
##'
##' # Postprocessing
##' plotGPareto(omEGO, add= FALSE, UQ_PF = TRUE, UQ_PS = TRUE, UQ_dens = TRUE)
##' 
##' #------------------------------------------------------------
##' # Example 2 : Surrogate-based multi-objective Optimization including a cheap function
##' #------------------------------------------------------------
##' set.seed(42)
##' library(DiceDesign)
##' 
##' d <- 2 
##' 
##' fname <- P1
##' n.grid <- 19
##' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
##' nappr <- 15 
##' design.grid <- maximinESE_LHS(lhsDesign(nappr, d, seed = 42)$design)$design
##' response.grid <- t(apply(design.grid, 1, fname))
##' 
##' mf1 <- km(~., design = design.grid, response = response.grid[,1])
##' mf2 <- km(~., design = design.grid, response = response.grid[,2])
##' model <- list(mf1, mf2)
##' 
##' nsteps <- 1 
##' lower <- rep(0, d)
##' upper <- rep(1, d)
##' 
##' # Optimization with fastfun: hypervolume with discrete search
##' 
##' optimcontrol <- list(method = "discrete", candidate.points = test.grid)
##' omEGO3 <- GParetoptim(model = model, fn = fname, cheapfn = branin, crit = "SMS",
##'                       nsteps = nsteps, lower = lower, upper = upper,
##'                       optimcontrol = optimcontrol)
##' print(omEGO3$par)
##' print(omEGO3$values) 
##' plotGPareto(omEGO3)
##' 
##' #------------------------------------------------------------
##' # Example 3 : Surrogate-based multi-objective Optimization (4 objectives)
##' #------------------------------------------------------------
##' set.seed(42)
##' library(DiceDesign)
##' 
##' d <- 5 
##' 
##' fname <- DTLZ3
##' nappr <- 25
##' design.grid <- maximinESE_LHS(lhsDesign(nappr, d, seed = 42)$design)$design
##' response.grid <- t(apply(design.grid, 1, fname, nobj = 4))
##' mf1 <- km(~., design = design.grid, response = response.grid[,1])
##' mf2 <- km(~., design = design.grid, response = response.grid[,2])
##' mf3 <- km(~., design = design.grid, response = response.grid[,3])
##' mf4 <- km(~., design = design.grid, response = response.grid[,4])
##' 
##' # Optimization
##' nsteps <- 5 
##' lower <- rep(0, d) 
##' upper <- rep(1, d)     
##' omEGO4 <- GParetoptim(model = list(mf1, mf2, mf3, mf4), fn = fname, crit = "EMI",
##'                       nsteps = nsteps, lower = lower, upper = upper, nobj = 4)
##' print(omEGO4$par)
##' print(omEGO4$values)
##' plotGPareto(omEGO4)
##' 
##' #------------------------------------------------------------
##' # Example 4 : quantification of uncertainty on Pareto front
##' #------------------------------------------------------------
##' library(DiceDesign)
##' set.seed(42)
##' 
##' nvar <- 2
##' 
##' # Test function P1
##' fname <- "P1"
##' 
##' # Initial design
##' nappr <- 10
##' design.grid <- maximinESE_LHS(lhsDesign(nappr, nvar, seed = 42)$design)$design
##' response.grid <- t(apply(design.grid, 1, fname))
##' 
##' PF <- t(nondominated_points(t(response.grid)))
##' 
##' # kriging models : matern5_2 covariance structure, linear trend, no nugget effect
##' mf1 <- km(~., design = design.grid, response = response.grid[,1])
##' mf2 <- km(~., design = design.grid, response = response.grid[,2])
##' 
##' # Conditional simulations generation with random sampling points 
##' nsim <- 100 # increase for better results
##' npointssim <- 1000 # increase for better results
##' Simu_f1 <- matrix(0, nrow = nsim, ncol = npointssim)
##' Simu_f2 <- matrix(0, nrow = nsim, ncol = npointssim)
##' design.sim <- array(0, dim = c(npointssim, nvar, nsim))
##' 
##' for(i in 1:nsim){
##'   design.sim[,,i] <- matrix(runif(nvar*npointssim), nrow = npointssim, ncol = nvar)
##'   Simu_f1[i,] <- simulate(mf1, nsim = 1, newdata = design.sim[,,i], cond = TRUE,
##'                          checkNames = FALSE, nugget.sim = 10^-8)
##'   Simu_f2[i,] <- simulate(mf2, nsim = 1, newdata = design.sim[,,i], cond = TRUE,
##'                          checkNames = FALSE, nugget.sim = 10^-8)
##' }
##' 
##' # Computation of the attainment function and Vorob'ev Expectation
##' CPF1 <- CPF(Simu_f1, Simu_f2, response.grid, paretoFront = PF)
##' 
##' summary(CPF1)
##' 
##' plot(CPF1)
##' 
##' # Display of the symmetric deviation function
##' plotSymDevFun(CPF1)
##' }
NULL 