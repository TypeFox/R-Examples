##' Find the design that maximizes the probability of dominating a target given by the user.
##'@title Get design corresponding to an objective target
##'@param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
##'@param target vector corresponding to the desired output in the objective space,
##'@param lower vector of lower bounds for the variables to be optimized over,
##'@param upper vector of upper bounds for the variables to be optimized over,
##'@param optimcontrol optional list of control parameters for optimization of the selected infill criterion. 
##'       "\code{method}" set the optimization method; one can 
##'        choose between "\code{discrete}", "\code{pso}" and "\code{genoud}". For each method, further parameters can be set.\cr 
##'        For "\code{discrete}", one has to provide the argument "\code{candidate.points}". \cr
##'        For "\code{pso}", one can control the maximum number of iterations "\code{maxit}" (\code{400}) and the population size "\code{s}"
##'        (default :  \code{max(20, floor(10+2*sqrt(length(dim))))} (see \code{\link[pso]{psoptim}}). \cr
##'        For "\code{genoud}", one can control, among others, "\code{pop.size}" (default :  \code{[N = 3*2^dim} for \code{dim < 6} and  \code{N = 32*dim} otherwise]),
##' "\code{max.generations}" (\code{12}), "\code{wait.generations}" (\code{2}), "\code{BFGSburnin}" (\code{2}), \code{BFGSmaxit} (\code{N}) and \code{solution.tolerance} (\code{1e-21})
##'  of function "\code{genoud}" (see \code{\link[rgenoud]{genoud}}). Numbers into brackets are the default values.

##'@return 
##' A list with components:
##' \itemize{
##' \item \code{par}: best design found,
##' \item \code{value}: probabilitity that the design dominates the target,
##' \item \code{mean}: kriging mean of the objectives at the design,
##' \item \code{sd}: prediction standard deviation at the design.
##' }
##'@export
##'@examples
##' \dontrun{
##' 
##' #---------------------------------------------------------------------------
##' # Example of interactive optimization
##' #---------------------------------------------------------------------------
##' 
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' d <- 2
##' n.obj <- 2 
##' fun <- "P1" 
##' n.grid <- 51
##' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
##' nappr <- 20 
##' design.grid <- round(maximinESE_LHS(lhsDesign(nappr, d, seed = 42)$design)$design, 1)
##' response.grid <- t(apply(design.grid, 1, fun))
##' paretoFront <- t(nondominated_points(t(response.grid)))
##' mf1 <- km(~., design = design.grid, response = response.grid[,1])
##' mf2 <- km(~., design = design.grid, response = response.grid[,2])
##' model <- list(mf1, mf2)
##' lower <- rep(0, d); upper <- rep(1, d)
##' 
##' sol <- GParetoptim(model, fun, crit = "SUR", nsteps = 5, lower = lower, upper = upper) 
##' 
##' plotGPareto(sol)
##' 
##' target1 <- c(15, -25)
##' points(x = target1[1], y = target1[2], col = "black", pch = 13)
##' 
##' nDesign <- getDesign(sol$lastmodel, target = target1, lower = rep(0, d), upper = rep(1, d))
##' points(t(nDesign$mean), col = "green", pch = 20)
##' 
##' target2 <- c(48, -27)
##' points(x = target2[1], y = target2[2], col = "black", pch = 13)
##' nDesign2 <- getDesign(sol$lastmodel, target = target2, lower = rep(0, d), upper = rep(1, d))
##' points(t(nDesign2$mean), col = "darkgreen", pch = 20)
##' }
getDesign <- function(model, target, lower, upper, optimcontrol = NULL){
  if(is.null(optimcontrol$method)) optimcontrol$method <- "pso"
  
  d <- model[[1]]@d
  
  if(optimcontrol$method=="discrete"){
    opt <- apply(optimcontrol$candidate.points, 1, prob.of.dominating, model = model, target = target)
    best <- which.max(opt)
    par <- optimcontrol$candidate.points[best,]
    value <- opt[best]
  }
  
  if(optimcontrol$method=="pso"){
    control <- list(fnscale=-1, maxit=optimcontrol$maxit, s = optimcontrol$s)
    if (is.null(control$maxit))   control$maxit=400
    if (is.null(control$s)) control$s = max(floor(10+2*sqrt(d)), 20)
    
    o <- psoptim(par = rep(NA, d), prob.of.dominating, lower = lower, upper = upper, control = control, model=model,
                 target = target)
    par <- t(as.matrix(o$par))
    colnames(par) <- colnames(model[[1]]@X)
    value <- as.matrix(-o$value)
    
  }  
  
  if(optimcontrol$method=="genoud"){
    if (d <= 6) N <- 3 * 2^d else N <- 32 * d
    
    if (is.null(optimcontrol$pop.size))         optimcontrol$pop.size <- N
    if (is.null(optimcontrol$max.generations))  optimcontrol$max.generations <- 12
    if (is.null(optimcontrol$wait.generations)) optimcontrol$wait.generations <- 2
    if (is.null(optimcontrol$BFGSburnin))       optimcontrol$BFGSburnin <- 2
    if (is.null(optimcontrol$parinit))          optimcontrol$parinit <- lower + runif(d) * (upper - lower)
    if (is.null(optimcontrol$unif.seed))        optimcontrol$unif.seed <- 1
    if (is.null(optimcontrol$int.seed))         optimcontrol$int.seed <- 1
    if (is.null(optimcontrol$print.level))           optimcontrol$print.level <- 0
    if (is.null(optimcontrol$BFGSmaxit))             optimcontrol$BFGSmaxit <- N
    if (is.null(optimcontrol$solution.tolerance))    optimcontrol$solution.tolerance <- 1e-21
    
    # Mutations
    if (is.null(optimcontrol$P1)) optimcontrol$P1<-50
    if (is.null(optimcontrol$P2)) optimcontrol$P2<-50
    if (is.null(optimcontrol$P3)) optimcontrol$P3<-50
    if (is.null(optimcontrol$P4)) optimcontrol$P4<-50
    if (is.null(optimcontrol$P5)) optimcontrol$P5<-50
    if (is.null(optimcontrol$P6)) optimcontrol$P6<-50
    if (is.null(optimcontrol$P7)) optimcontrol$P7<-50
    if (is.null(optimcontrol$P8)) optimcontrol$P8<-50
    if (is.null(optimcontrol$P9)) optimcontrol$P9<-0
    
    domaine <- cbind(lower, upper)
    
    o <- genoud(fn=prob.of.dominating, nvars=d, max=TRUE, pop.size=optimcontrol$pop.size,
                max.generations=optimcontrol$max.generations,wait.generations=optimcontrol$wait.generations,
                hard.generation.limit=TRUE, starting.values=optimcontrol$parinit, MemoryMatrix=TRUE,
                Domains=domaine, default.domains=10, solution.tolerance=optimcontrol$solution.tolerance,
                boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE, BFGS=TRUE,
                data.type.int=FALSE, hessian=FALSE, unif.seed=optimcontrol$unif.seed, 
                int.seed=optimcontrol$int.seed,print.level=optimcontrol$print.level, share.type=0, instance.number=0,
                output.path="stdout", output.append=FALSE, project.path=NULL,
                P1=optimcontrol$P1, P2=optimcontrol$P2, P3=optimcontrol$P3, 
                P4=optimcontrol$P4, P5=optimcontrol$P5, P6=optimcontrol$P6,
                P7=optimcontrol$P7, P8=optimcontrol$P8, P9=optimcontrol$P9,
                P9mix=NULL, BFGSburnin=optimcontrol$BFGSburnin,BFGSfn=NULL, BFGShelp=NULL,
                control = list(maxit = optimcontrol$BFGSmaxit),
                cluster=FALSE, balance=FALSE, debug=FALSE,
                model=model, target = target)
    
    par <- t(as.matrix(o$par))
    colnames(par) <- colnames(model[[1]]@X)
    value <- as.matrix(o$value)
  } 
  
  # Prediction at this best design
  tmp <- lapply(model, predict, newdata = par, type = "UK",
                light.return = TRUE, checkNames = F)
  
  tmp <- matrix(unlist(tmp), ncol=5, byrow=TRUE)
  
  return(list(par=par, value=value, mean = tmp[,2], sd = tmp[,3])) #,
              #lower95 = tmp[,4], upper95 = tmp[,5]))
}

prob.of.dominating <- function(X, model, target){
  tmp <- lapply(model, predict, newdata = data.frame(t(X)), type = "UK",
                light.return = TRUE, checkNames = F)
  prob <- prod(mapply(univariateProb, list = tmp, target = target))
  return(prob)
}

univariateProb <- function(list, target){
  return(pnorm(q = target, mean = list$mean, sd = list$sd))
}