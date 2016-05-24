##' Computes a slightly modified infill Criterion of the SMS-EGO.
##' To avoid numerical instabilities, an additional penalty is added to the new point if it is too close to an existing observation.
##' @title Analytical expression of the SMS-EGO criterion with m>1 objectives
##' @param x a vector representing the input for which one wishes to calculate the criterion,
##' @param model a list of objects of class \code{\link[DiceKriging]{km}} (one for each objective),
##' @param paretoFront (optional) matrix corresponding to the Pareto front of size \code{[n.pareto x n.obj]}, 
##' @param critcontrol list with arguments: 
##' \itemize{
##' \item \code{currentHV} current hypervolume;
##' \item \code{refPoint} reference point for hypervolume computations;
##' \item \code{epsilon} optional value to use in additive epsilon dominance;
##' \item \code{gain} optional gain factor for sigma.
##' }
##'        Options for the \code{\link[GPareto]{checkPredict}} function: \code{threshold} (\code{1e-4}) and \code{distance} (\code{covdist}) are used to avoid numerical issues occuring when adding points too close to the existing ones.
##' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend 
##'        estimation has to be taken into account.
##' @references 
##' W. Ponweiser, T. Wagner, D. Biermann, M. Vincze (2008), Multiobjective Optimization on a Limited Budget of Evaluations Using Model-Assisted S-Metric Selection,
##' \emph{Parallel Problem Solving from Nature}, pp. 784-794. Springer, Berlin. \cr \cr
##' T. Wagner, M. Emmerich, A. Deutz, W. Ponweiser (2010), On expected-improvement criteria for model-based multi-objective optimization.   
##' \emph{Parallel Problem Solving from Nature}, pp. 718-727. Springer, Berlin.
##' @seealso \code{\link[GPareto]{crit_EHI}}, \code{\link[GPareto]{crit_SUR}}, \code{\link[GPareto]{crit_EMI}}.
##' @return Value of the criterion.
##' @examples
##' #---------------------------------------------------------------------------
##' # SMS-EGO surface associated with the "P1" problem at a 15 points design
##' #---------------------------------------------------------------------------
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' n_var <- 2 
##' f_name <- "P1" 
##' n.grid <- 26
##' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
##' n_appr <- 15 
##' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
##' response.grid <- t(apply(design.grid, 1, f_name))
##' PF <- t(nondominated_points(t(response.grid)))
##' mf1 <- km(~., design = design.grid, response = response.grid[,1])
##' mf2 <- km(~., design = design.grid, response = response.grid[,2])
##' 
##' model <- list(mf1, mf2)
##' critcontrol <- list(refPoint = c(300, 0), currentHV = dominated_hypervolume(t(PF), c(300, 0)))
##' SMSEGO_grid <- apply(test.grid, 1, crit_SMS, model = model,
##'                      paretoFront = PF, critcontrol = critcontrol)
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid),
##'                matrix(pmax(0, SMSEGO_grid), nrow = n.grid), nlevels = 50,
##'                main = "SMS-EGO criterion (positive part)", xlab = expression(x[1]),
##'                ylab = expression(x[2]), color = terrain.colors,
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1],design.grid[,2], pch = 21, bg = "white")
##'                             }
##'               )
##' @export

crit_SMS <- function(x, model, paretoFront=NULL, critcontrol=list(epsilon=1e-6, gain=1), type="UK")
{
  # Slightly modified infill Criterion of the SMSEGO
  # Ponweiser, Wagner et al. (Proc. 2008 PPSN, pp. 784-794)
  # ***********************************************************************
  # arguments
  # x:             decision vector to be evaluated
  # model:         d-dimensional cell array of models for each objective
  # ref:           d-dimensional anti-ideal reference point
  # paretoFront:   current Pareto front approximation
  # currentHV:     hypervolume of current front with respect to ref
  # epsilon:       epsilon to use in additive epsilon dominance
  # gain:          gain factor for sigma (optional)
  
  n.obj <- length(model)
  d <- model[[1]]@d
  x.new <- matrix(x, 1, d)
  
  distp <- 0  # penalty if too close in the checkPredict sense
  
  if(checkPredict(x, model, type = type, distance = critcontrol$distance, threshold = critcontrol$threshold)){
    # return(0) Not compatible with penalty with SMS
    distp <- 1 # may be changed
  }#else{
  
  refPoint  <- critcontrol$refPoint
  currentHV <- critcontrol$currentHV
  epsilon   <- critcontrol$epsilon
  gain      <- critcontrol$gain
  
  if(is.null(paretoFront) || is.null(refPoint)) {
    observations <- matrix(0, model[[1]]@n, n.obj)
    for (i in 1:n.obj) observations[,i] <- model[[i]]@y
  }
  if(is.null(paretoFront)) paretoFront <- t(nondominated_points(t(observations)))
  if (is.null(refPoint)){
    refPoint    <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj)
    cat("No refPoint provided, ", signif(refPoint, 3), "used \n")
  }
  
  if (is.null(currentHV)) currentHV <- dominated_hypervolume(points=t(paretoFront), ref=refPoint)
  if (is.null(epsilon))   epsilon <- 1e-6
  if (is.null(gain))      gain <- 1
  n.pareto <- nrow(paretoFront)
  
  mu    <- rep(NaN, n.obj)
  sigma <- rep(NaN, n.obj)
  for (i in 1:n.obj){    
    pred     <- predict(object=model[[i]], newdata=x.new, type=type, checkNames = FALSE)
    mu[i]    <- pred$mean
    sigma[i] <- pred$sd
  }
  potSol <- mu - gain*sigma
  penalty <- distp
  for (j in 1:n.pareto){
    # assign penalty to all epsilon-dominated solutions
    if (min(paretoFront[j,] <= potSol + epsilon)){
      p <- -1 + prod(1 + pmax(potSol - paretoFront[j,], rep(0, n.obj)))
      penalty <- max(penalty, p)
    }
  }
  if (penalty == 0){
    # non epsilon-dominated solution
    potFront <- rbind(paretoFront, potSol)
    mypoints <- t(nondominated_points(t(potFront)))
    if (is.null(nrow(mypoints))) mypoints <- matrix(mypoints, 1, n.obj)
    myhv <- dominated_hypervolume(points=t(mypoints), ref=refPoint)
    f    <- currentHV - myhv
  } else{
    f <- penalty
  }
  #}
  
  return(-f)
}