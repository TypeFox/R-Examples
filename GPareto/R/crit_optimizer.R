##' Given a list of objects of class \code{\link[DiceKriging]{km}} and a set of tuning
##' parameters (\code{lower, upper and critcontrol}), \code{crit_optimizer} performs
##' the maximization of an Expected Improvement or SUR criterion and delivers
##' the next point to be visited in a multi-objective EGO-like procedure. \cr \cr
##' The latter maximization relies either on a genetic algorithm using derivatives,
##' \code{\link[rgenoud]{genoud}}, particle swarm algorithm \code{\link[pso]{pso}},
##'  exhaustive search at pre-specified points or on a user defined method. 
##' It is important to remark that the information
##' needed about the objective function reduces here to the vector of response values
##' embedded in the models (no call to the objective functions or simulators (except with \code{cheapfn})).
##'
##' @title Maximization of multiobjective Expected Improvement criteria
##' 
##' @param crit sampling criterion. Four choices are available : "\code{SMS}", "\code{EHI}", "\code{EMI}" and "\code{SUR}",
##' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
##' @param lower vector of lower bounds for the variables to be optimized over,
##' @param upper vector of upper bounds for the variables to be optimized over,
##' @param cheapfn optional additional fast-to-evaluate objective function (handled next with class \code{\link[GPareto]{fastfun}}), which does not need a kriging model,
##' @param type "\code{SK}" or "\code{UK}" (default), depending whether uncertainty related to trend estimation has to be taken into account.
##' @param paretoFront (optional) matrix corresponding to the Pareto front of size \code{[n.pareto x n.obj]}, 
##' @param critcontrol optional list of control parameters for criterion \code{crit}, see details.
##' Options for the \code{\link[GPareto]{checkPredict}} function: \code{threshold} (\code{1e-4}) and \code{distance} (\code{covdist}) are used to avoid numerical issues occuring when adding points too close to the existing ones.
##' @param optimcontrol optional list of control parameters for optimization of the selected infill criterion. 
##'       "\code{method}" set the optimization method; one can 
##'        choose between "\code{discrete}", "\code{pso}" and "\code{genoud}" or a user defined method name (passed to \code{\link[base]{match.fun}}). For each method, further parameters can be set.\cr 
##'        For "\code{discrete}", one has to provide the argument "\code{candidate.points}". \cr
##'        For "\code{pso}", one can control the maximum number of iterations "\code{maxit}" (\code{400}) and the population size "\code{s}"
##'        (default :  \code{max(20, floor(10+2*sqrt(length(dim))))} (see \code{\link[pso]{psoptim}}). \cr
##'        For "\code{genoud}", one can control, among others, "\code{pop.size}" (default :  \code{[N = 3*2^dim} for \code{dim < 6} and  \code{N = 32*dim} otherwise]),
##' "\code{max.generations}" (\code{12}), "\code{wait.generations}" (\code{2}), "\code{BFGSburnin}" (\code{2}), \code{BFGSmaxit} (\code{N}) and \code{solution.tolerance} (\code{1e-21})
##'  of function "\code{genoud}" (see \code{\link[rgenoud]{genoud}}). Numbers into brackets are the default values.\cr
##'  For a user defined method, it must have arguments like the default \code{\link[stats]{optim}} method, i.e. \code{par}, \code{fn}, \code{lower}, \code{upper}, \code{...} and eventually \code{control}, and return a list with \code{par} and \code{value}.
##' @param nsteps.remaining Number of iterations remaining in the optimization loop. 
##'        Used by "\code{SMS}" to determine the parameter epsilon.
##'        
##' @return A list with components:  
##'  \itemize{
##'  \item{\code{par}}{: The best set of parameters found,}
##'  \item{\code{value}}{: The value of expected improvement at \code{par}.}
##'  }
##' 
##' @details
##' Extension of the function \code{\link[DiceOptim]{max_EI}} for multi-objective optimization.\cr
##' Available infill criteria with \code{crit} are : \cr
##' \itemize{
##' \item Expected Hypervolume Improvement (\code{EHI}) \code{\link[GPareto]{crit_EHI}},
##' \item SMS criterion (\code{SMS}) \code{\link[GPareto]{crit_SMS}},
##' \item Expected Maximin Improvement (\code{EMI}) \code{\link[GPareto]{crit_EMI}},
##' \item Stepwise Uncertainty Reduction of the excursion volume (\code{SUR}) \code{\link[GPareto]{crit_SUR}}
##' }
##' Depending on the selected criterion, parameters such as a reference point for \code{SMS} and \code{EHI} or arguments for \code{\link[GPareto]{integration_design_optim}} with \code{SUR} can be given with \code{critcontrol}.
##' Also options for \code{\link[GPareto]{checkPredict}} are available.
##' More precisions are given in the corresponding help pages. 
##' @importFrom pso psoptim
##' @importFrom KrigInv precomputeUpdateData
##' @references
##' W.R. Jr. Mebane and J.S. Sekhon (2011), Genetic optimization using derivatives: The rgenoud package for R, \emph{Journal of Statistical Software}, 42(11), 1-26 \cr
##' 
##' D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global optimization of expensive black-box functions, \emph{Journal of Global Optimization}, 13, 455-492.
##' @examples
##' \dontrun{
##' #---------------------------------------------------------------------------
##' # EHI surface associated with the "P1" problem at a 15 points design
##' #---------------------------------------------------------------------------
##' 
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' d <- 2
##' n.obj <- 2 
##' fname <- "P1" 
##' n.grid <- 51
##' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
##' nappr <- 15 
##' design.grid <- round(maximinESE_LHS(lhsDesign(nappr, d, seed = 42)$design)$design, 1)
##' response.grid <- t(apply(design.grid, 1, fname))
##' paretoFront <- t(nondominated_points(t(response.grid)))
##' mf1 <- km(~., design = design.grid, response = response.grid[,1])
##' mf2 <- km(~., design = design.grid, response = response.grid[,2])
##' model <- list(mf1, mf2)
##' 
##' EHI_grid <- apply(test.grid, 1, crit_EHI, model = list(mf1, mf2),
##'                   critcontrol = list(refPoint = c(300, 0)))
##'
##' lower <- rep(0, d)
##' upper <- rep(1, d)     
##' 
##' omEGO <- crit_optimizer(crit = "EHI", model = model,  lower = lower, upper = upper, 
##'                  optimcontrol = list(method = "genoud", pop.size = 200, BFGSburnin = 2),
##'                  critcontrol = list(refPoint = c(300, 0)))
##'                  
##' print(omEGO)
##'  
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(EHI_grid, nrow = n.grid), main = "Expected Hypervolume Improvement", 
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors,
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[, 1], design.grid[, 2], pch = 21, bg = "white");
##'                             points(omEGO$par, col = "red", pch = 4)
##'                             }
##'               )
##'               
##' #---------------------------------------------------------------------------
##' # SMS surface associated with the "P1" problem at a 15 points design
##' #---------------------------------------------------------------------------
##' 
##' 
##' SMS_grid <- apply(test.grid, 1, crit_SMS, model = model,
##'                   critcontrol = list(refPoint = c(300, 0)))
##'
##' lower <- rep(0, d)
##' upper <- rep(1, d)     
##' 
##' omEGO2 <- crit_optimizer(crit = "SMS", model = model,  lower = lower, upper = upper, 
##'                  optimcontrol = list(method="genoud", pop.size = 200, BFGSburnin = 2),
##'                  critcontrol = list(refPoint = c(300, 0)))
##'                  
##' print(omEGO2)
##'  
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(pmax(0,SMS_grid), nrow = n.grid), main = "SMS Criterion (>0)",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors,
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[, 1], design.grid[, 2], pch = 21, bg = "white");
##'                             points(omEGO2$par, col = "red", pch = 4)
##'                             }
##'               )
##' #---------------------------------------------------------------------------
##' # Maximin Improvement surface associated with the "P1" problem at a 15 points design
##' #---------------------------------------------------------------------------
##' 
##' 
##' EMI_grid <- apply(test.grid, 1, crit_EMI, model = model,
##'                   critcontrol = list(nb_samp = 20, type ="EMI"))
##'
##' lower <- rep(0, d)
##' upper <- rep(1, d)     
##' 
##' omEGO3 <- crit_optimizer(crit = "EMI", model = model,  lower = lower, upper = upper, 
##'                  optimcontrol = list(method = "genoud", pop.size = 200, BFGSburnin = 2))
##'                  
##' print(omEGO3)
##'  
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(EMI_grid, nrow = n.grid), main = "Expected Maximin Improvement",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
##'                plot.axes = {axis(1);axis(2);
##'                             points(design.grid[, 1], design.grid[, 2], pch = 21, bg = "white");
##'                             points(omEGO3$par, col = "red", pch = 4)
##'                             }
##'               )
##' #---------------------------------------------------------------------------
##' # crit_SUR surface associated with the "P1" problem at a 15 points design
##' #---------------------------------------------------------------------------
##' library(KrigInv)
##' 
##' integration.param <- integration_design_optim(lower = c(0, 0), upper = c(1, 1), model = model)
##' integration.points <- as.matrix(integration.param$integration.points)
##' integration.weights <- integration.param$integration.weights
##' 
##' precalc.data <- list()
##' mn.X <- sn.X <- matrix(0, n.obj, nrow(integration.points))
##' 
##' for (i in 1:n.obj){
##'   p.tst.all <- predict(model[[i]], newdata = integration.points, type = "UK",
##'                        checkNames = FALSE)
##'   mn.X[i,] <- p.tst.all$mean
##'   sn.X[i,]   <- p.tst.all$sd
##'   precalc.data[[i]] <- precomputeUpdateData(model[[i]], integration.points)
##' }
##' critcontrol <- list(integration.points = integration.points,
##'                     integration.weights = integration.weights,
##'                     mn.X = mn.X, sn.X = sn.X, precalc.data = precalc.data)
##' EEV_grid <- apply(test.grid, 1, crit_SUR, model=model, paretoFront = paretoFront,
##'                   critcontrol = critcontrol)
##' lower <- rep(0, d)
##' upper <- rep(1, d)
##'                   
##' omEGO4 <- crit_optimizer(crit = "SUR", model = model,  lower = lower, upper = upper, 
##'                  optimcontrol = list(method = "genoud", pop.size = 200, BFGSburnin = 2))
##' print(omEGO4)
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid),
##'                matrix(pmax(0,EEV_grid), n.grid), main = "EEV criterion", nlevels = 50,
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors,
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             points(omEGO4$par, col = "red", pch = 4)
##'                             }
##'               )
##' 
##' # example using user defined optimizer, here L-BFGS-B from base optim
##' userOptim <- function(par, fn, lower, upper, control, ...){
##'   return(optim(par = par, fn = fn, method = "L-BFGS-B", lower = lower, upper = upper,
##'          control = control, ...))
##' }
##' omEGO4bis <- crit_optimizer(crit = "SUR", model = model,  lower = lower, upper = upper, 
##'                  optimcontrol = list(method = "userOptim"))
##' print(omEGO4bis)
##' 
##' 
##' #---------------------------------------------------------------------------
##' # crit_SMS surface with problem "P1" with 15 design points, using cheapfn
##' #---------------------------------------------------------------------------
##' 
##' # Optimization with fastfun: SMS with discrete search
##' # Separation of the problem P1 in two objectives: 
##' # the first one to be kriged, the second one with fastobj
##' 
##' # Definition of the fastfun
##' f2 <-   function(x){
##'   return(P1(x)[2])
##' }
##'  
##' SMS_grid_cheap <- apply(test.grid, 1, crit_SMS, model = list(mf1, fastfun(f2, design.grid)),
##'                         paretoFront = paretoFront, critcontrol = list(refPoint = c(300, 0)))
##' 
##' 
##' optimcontrol <- list(method = "pso")
##' model2 <- list(mf1)
##' omEGO5 <- crit_optimizer(crit = "SMS", model = model2,  lower = lower, upper = upper,
##'                          cheapfn = f2, critcontrol = list(refPoint = c(300, 0)),
##'                          optimcontrol = list(method = "genoud", pop.size = 200, BFGSburnin = 2))
##' print(omEGO5)
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                matrix(pmax(0, SMS_grid_cheap), nrow = n.grid), nlevels = 50,
##'                main = "SMS criterion with cheap 2nd objective (>0)", xlab = expression(x[1]),
##'                ylab = expression(x[2]), color = terrain.colors,
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             points(omEGO5$par, col = "red", pch = 4)
##'                             }
##'               )
##' }
##' @export

crit_optimizer <- function(crit = "SMS", model, lower, upper, cheapfn = NULL, type = "UK", paretoFront = NULL, 
                        critcontrol = NULL, optimcontrol = NULL, nsteps.remaining = 1){
  ###########################################################################################
  # Finds the maximizer of the criterion
  ###########################################################################################
  if(is.null(optimcontrol$method)) optimcontrol$method <- "genoud"
  
  if (!is.null(cheapfn)) {
    fastobs <- apply(model[[1]]@X, 1, cheapfn)
    fastmod <- fastfun(fn = cheapfn, design = model[[1]]@X, response = fastobs)
    model[[length(model)+1]] <- fastmod
  }
  
  d     <- model[[1]]@d
  n.obj <- length(model)

  
  #if (is.null(observations)) {
  observations <- matrix(0, model[[1]]@n, n.obj)
  for (i in 1:n.obj) observations[,i] <- model[[i]]@y
  #}
  
  if (is.null(paretoFront)) paretoFront <- t(nondominated_points(t(observations)))
  n.pareto <- nrow(paretoFront)
  if(is.null(critcontrol)) 
    critcontrol <- list()
  
  
  if (is.null(critcontrol$refPoint)){    
    critcontrol$refPoint <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj) ### May be changed
    if(crit == "SMS" | crit =="EHI")
      cat("No refPoint provided, ", signif(critcontrol$refPoint, 3), "used \n")
  } 
  # Different calls for crit_optimizer, depending on the criterion chosen
  if (crit=="SMS"){
    #-------------------------------------------------------
    currentHV <- dominated_hypervolume(points=t(paretoFront), ref=critcontrol$refPoint)
    if (n.pareto < 2){
      epsilon <- rep(0, n.obj)
    } else {
      spread <- apply(paretoFront,2,max) - apply(paretoFront,2,min)
      c <- 1 - (1 / (2^n.obj) )
      epsilon <- spread / (n.pareto + c * (nsteps.remaining-1))
    }
    gain <- -qnorm( 0.5*(0.5^(1/n.obj)) )
    criterion <- crit_SMS
    critcontrol   <- c(critcontrol, list(gain=gain, epsilon=epsilon, currentHV=currentHV))
    #-------------------------------------------------------
  } else if (crit=="EHI"){
    criterion <- crit_EHI
    
    if (is.unsorted(paretoFront[,1])){
      paretoFront <- paretoFront[sort(paretoFront[,1], index.return=TRUE)[[2]],,drop=FALSE]
    }
    if(n.obj > 2)
      critcontrol$seed <- sample.int(1e9,1)
    
    #-------------------------------------------------------
  } else if (crit == "EMI"){
    criterion <- crit_EMI
    critcontrol$seed <- sample.int(1e9,1)
    
    #-------------------------------------------------------
  } else if (crit=="SUR"){
    
    ## if values for critcontrol are not provided correctly (integration.points,mn.X,...)
    ## -> integration_design_optim is applied
    
    if(is.null(critcontrol$integration.points)){
      integration.param   <- integration_design_optim(critcontrol$SURcontrol, d, lower, upper, model=model)
      integration.points  <- as.matrix(integration.param$integration.points)
      integration.weights  <- integration.param$integration.weights
      
      precalc.data <- vector("list", n.obj)
      intpoints.oldmean <- intpoints.oldsd <- matrix(0, n.obj, nrow(integration.points))
      
      for (i in 1:n.obj){
        p.tst <- predict(model[[i]], newdata=integration.points, type=type, checkNames=FALSE)
        intpoints.oldmean[i,] <- p.tst$mean
        intpoints.oldsd[i,]   <- p.tst$sd
        if (max(intpoints.oldsd[i,]) != 0) precalc.data[[i]]     <- precomputeUpdateData(model[[i]], integration.points)
      }
      
      critcontrol <- c(critcontrol, list(integration.points=integration.points, integration.weights=integration.param$integration.weights,
                                         mn.X=intpoints.oldmean, sn.X=intpoints.oldsd, precalc.data=precalc.data))
    }

    
    ## Reorder current Pareto front if needed
    if (n.obj==2){ if (is.unsorted(paretoFront[,1]))  paretoFront <- paretoFront[sort(paretoFront[,1], index.return=TRUE)[[2]],,drop=FALSE]
    } else {       if (is.unsorted(-paretoFront[,3])) paretoFront <- paretoFront[sort(paretoFront[,3], index.return=TRUE, decreasing=TRUE)[[2]],,drop=FALSE]
    }

    criterion <- crit_SUR

    #-------------------------------------------------------
  }
  
  ########################################################################################
  ## Discrete Optimisation
  ########################################################################################
  
  if(optimcontrol$method=="discrete"){
    optim.points <- optimcontrol$candidate.points
    colnames(optim.points) <- colnames(model[[1]]@X)
    n.optim.points <- nrow(optim.points)
    all.crit <- rep(0, n.optim.points)
    critcontroldiscrete <- critcontrol
    
    for (k in 1:n.optim.points){
      if (crit=="SUR"){
        critcontroldiscrete$integration.points <- critcontrol$integration.points[-k,,drop=FALSE]
        if (!is.null(critcontroldiscrete$integration.weights)) critcontroldiscrete$integration.weights <- critcontrol$integration.weights[-k]
        critcontroldiscrete$mn.X <- critcontrol$mn.X[,-k,drop=FALSE]
        critcontroldiscrete$sn.X <- critcontrol$sn.X[,-k,drop=FALSE]
        for (i in 1:n.obj){
          if (!is.null(precalc.data[[i]])) {
          critcontroldiscrete$precalc.data[[i]]$Kinv.c.olddata <- critcontrol$precalc.data[[i]]$Kinv.c.olddata[,-k,drop=FALSE]
          critcontroldiscrete$precalc.data[[i]]$first.member   <- critcontrol$precalc.data[[i]]$first.member[-k]
          }
        }
      }
      all.crit[k] <- criterion(x=as.numeric(optim.points[k,]), paretoFront=paretoFront, model=model, type=type,
                               critcontrol=critcontroldiscrete)
    }
    value <- max(all.crit)
    i.best <- which(all.crit==value)[sample.int(length(which(all.crit==value)), 1)]
    par <- matrix(optim.points[i.best,],nrow=1, ncol=model[[1]]@d)
    colnames(par) <- colnames(model[[1]]@X)
    
    return(list(par=par, value=value, index=i.best))
  }
  
  ########################################################################################
  ## Optimization with Genoud
  ########################################################################################
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
    
    o <- genoud(fn=criterion, nvars=d, max=TRUE, pop.size=optimcontrol$pop.size,
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
                model=model, type=type, paretoFront=paretoFront, 
                critcontrol=critcontrol)
    
    par <- t(as.matrix(o$par))
    colnames(par) <- colnames(model[[1]]@X)
    value <- as.matrix(o$value)
    return(list(par=par, value=value))
  }
  ########################################################################################
  ## Optimization with PSO
  ########################################################################################
  if(optimcontrol$method=="pso"){
    
    control <- list(fnscale=-1, maxit=optimcontrol$maxit, s = optimcontrol$s)
    if (is.null(control$maxit))   control$maxit=400
    if (is.null(control$s)) control$s = max(floor(10+2*sqrt(d)), 20)
    
    o <- psoptim(par = rep(NA, d) , criterion, lower = lower, upper = upper, control = control, model=model, type=type, 
                 paretoFront=paretoFront, critcontrol=critcontrol)
    par <- t(as.matrix(o$par))
    colnames(par) <- colnames(model[[1]]@X)
    value <- as.matrix(o$value)
    return(list(par=par, value=value))
  }
  
  ########################################################################################
  ## Used defined optimization
  ########################################################################################
  optimMethodUser <- match.fun(optimcontrol$method)
  control <- optimcontrol
  control$method <- NULL # avoid warnings
  control$notrace <- NULL # avoid warnings
  rand_start <- runif(d)*(upper-lower) + lower # random start point
  o <- optimMethodUser(par = rand_start, fn = criterion, lower = lower, upper = upper,
                       control = control, model=model, type=type, 
                       paretoFront=paretoFront, critcontrol=critcontrol)
  par <- t(as.matrix(o$par))
  colnames(par) <- colnames(model[[1]]@X)
  value <- as.matrix(o$value)
  return(list(par=par, value=value))
}