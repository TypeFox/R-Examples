##' Executes \code{nsteps} iterations of multi-objective EGO methods to objects of class \code{\link[DiceKriging]{km}}.
##' At each step, kriging models are re-estimated (including covariance parameters re-estimation)
##'  based on the initial design points plus the points visited during all previous iterations;
##'  then a new point is obtained by maximizing one of the four multi-objective Expected Improvement criteria available. 
##' @title Sequential multi-objective Expected Improvement maximization and model re-estimation,
##'  with a number of iterations fixed in advance by the user
##' @details Extension of the function \code{\link[DiceOptim]{EGO.nsteps}} for multi-objective optimization.\cr
##' Available infill criteria with \code{crit} are: \cr
##' \itemize{
##' \item Expected Hypervolume Improvement (\code{EHI}) \code{\link[GPareto]{crit_EHI}},
##' \item SMS criterion (\code{SMS}) \code{\link[GPareto]{crit_SMS}},
##' \item Expected Maximin Improvement (\code{EMI}) \code{\link[GPareto]{crit_EMI}},
##' \item Stepwise Uncertainty Reduction of the excursion volume (\code{SUR}) \code{\link[GPareto]{crit_SUR}}.
##' }
##' Depending on the selected criterion, parameters such as reference point for \code{SMS} and \code{EHI} or arguments for \code{\link[GPareto]{integration_design_optim}} with \code{SUR} can be given with \code{critcontrol}.
##' Also options for \code{\link[GPareto]{checkPredict}} are available.
##' More precisions are given in the corresponding help pages.\cr
##' 
##' Display of results and various post-processings are available with \code{\link[GPareto]{plotGPareto}}.  
##' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
##' @param fn the multi-objective function to be minimized (vectorial output), found by a call to \code{\link[base]{match.fun}},
##' @param cheapfn optional additional fast-to-evaluate objective function (handled next with class \code{\link[GPareto]{fastfun}}), which does not need a kriging model, handled by a call to \code{\link[base]{match.fun}},
##' @param crit choice of multi-objective improvement function: "\code{SMS}", "\code{EHI}", "\code{EMI}" or "\code{SUR}",
##' see details below,
##' @param nsteps an integer representing the desired number of iterations,
##' @param lower vector of lower bounds for the variables to be optimized over,
##' @param upper vector of upper bounds for the variables to be optimized over,
##' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend estimation has to be taken into account,
##' @param cov.reestim optional boolean specifying if the kriging hyperparameters should be re-estimated at each iteration,
##' @param critcontrol optional list of parameters for criterion \code{crit}, see details,
##' @param optimcontrol an optional list of control parameters for optimization of the selected infill criterion:
##' "\code{method}" can be set to "\code{discrete}", "\code{pso}", "\code{genoud}" or a user defined method name (passed to \code{\link[base]{match.fun}}). For "\code{discrete}", a matrix \code{candidate.points} must be given.
##' For "\code{pso}" and "\code{genoud}", specific parameters to the chosen method can also be specified  (see \code{\link[rgenoud]{genoud}} and \code{\link[pso]{psoptim}}).
##' A user defined method must have arguments like the default \code{\link[stats]{optim}} method, i.e. \code{par}, \code{fn}, \code{lower}, \code{upper}, \code{...} and eventually \code{control}.\cr
##' Option \code{notrace} can be set to \code{TRUE} to suppress printing of the optimization progresses. 
##' 
##' @param ... additional parameters to be given to the objective \code{fn}.
##' @export
##' @return
##' A list with components:
##' \itemize{
##' \item{\code{par}}{: a data frame representing the additional points visited during the algorithm,}
##' \item{\code{values}}{: a data frame representing the response values at the points given in \code{par},}
##' \item{\code{nsteps}}{: an integer representing the desired number of iterations (given in argument),}
##' \item{\code{lastmodel}}{: a list of objects of class \code{\link[DiceKriging]{km}} corresponding to the last kriging models fitted.}
##' If a problem occurs during either model updates or criterion maximization, the last working model and corresponding values are returned.
##' }
##' 
##' @references 
##' M. T. Emmerich, A. H. Deutz, J. W. Klinkenberg (2011), Hypervolume-based expected improvement: Monotonicity properties and exact computation,
##' \emph{Evolutionary Computation (CEC)}, 2147-2154. \cr \cr
##' V. Picheny (2014), Multiobjective optimization using Gaussian process emulators via stepwise uncertainty reduction, 
##' \emph{Statistics and Computing}. \cr \cr
##' T. Wagner, M. Emmerich, A. Deutz, W. Ponweiser (2010), On expected-improvement criteria for model-based multi-objective optimization.   
##' \emph{Parallel Problem Solving from Nature}, 718-727, Springer, Berlin. \cr \cr
##' J. D. Svenson (2011), \emph{Computer Experiments: Multiobjective Optimization and Sensitivity Analysis}, Ohio State university, PhD thesis. 
##' 
##' @importFrom stats runif pnorm qnorm
##' 
##' @examples
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' d <- 2 
##' 
##' fname <- ZDT3
##' n.grid <- 21
##' test.grid = expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
##' nappr <- 15 
##' design.grid <- maximinESE_LHS(lhsDesign(nappr, d, seed = 42)$design)$design
##' response.grid <- t(apply(design.grid, 1, fname))
##' Front_Pareto <- t(nondominated_points(t(response.grid)))
##' 
##' mf1 <- km(~., design = design.grid, response = response.grid[, 1])
##' mf2 <- km(~., design = design.grid, response = response.grid[, 2])
##' model <- list(mf1, mf2)
##' 
##' nsteps <- 2
##' lower <- rep(0, d)
##' upper <- rep(1, d)
##' 
##' # Optimization 1: EHI with pso
##' optimcontrol <- list(method = "pso", maxit = 20)
##' critcontrol <- list(refPoint = c(1, 10))
##' omEGO1 <- GParetoptim(model = model, fn = fname, crit = "EHI", nsteps = nsteps,
##'                      lower = lower, upper = upper, critcontrol = critcontrol,
##'                      optimcontrol = optimcontrol)
##' print(omEGO1$par)
##' print(omEGO1$values)
##' 
##' \dontrun{
##' # Optimization 2: SMS with discrete search
##' optimcontrol <- list(method = "discrete", candidate.points = test.grid)
##' critcontrol <- list(refPoint = c(1, 10))
##' omEGO2 <- GParetoptim(model = model, fn = fname, crit = "SMS", nsteps = nsteps,
##'                      lower = lower, upper = upper, critcontrol = critcontrol,
##'                      optimcontrol = optimcontrol)
##' print(omEGO2$par)
##' print(omEGO2$values)
##' 
##' # Optimization 3: SUR with genoud
##' optimcontrol <- list(method = "genoud", pop.size = 20, max.generations = 10)
##' critcontrol <- list(SURcontrol = list(distrib = "SUR", n.points = 100))
##' omEGO3 <- GParetoptim(model = model, fn = fname, crit = "SUR", nsteps = nsteps,
##'                      lower = lower, upper = upper, critcontrol = critcontrol,
##'                      optimcontrol = optimcontrol)
##' print(omEGO3$par)
##' print(omEGO3$values)
##' 
##' # Optimization 4: EMI with pso
##' optimcontrol <- list(method = "pso", maxit = 20)
##' critcontrol <- list(nbsamp = 200)
##' omEGO4 <- GParetoptim(model = model, fn = fname, crit = "EMI", nsteps = nsteps,
##'                      lower = lower, upper = upper, optimcontrol = optimcontrol)
##' print(omEGO4$par)
##' print(omEGO4$values)
##' 
##' # graphics
##' sol.grid <- apply(expand.grid(seq(0, 1, length.out = 100),
##'                               seq(0, 1, length.out = 100)), 1, fname)
##' plot(t(sol.grid), pch = 20, col = rgb(0, 0, 0, 0.05), xlim = c(0, 1),
##'      ylim = c(-2, 10), xlab = expression(f[1]), ylab = expression(f[2]))
##' plotGPareto(list = omEGO1, add = TRUE,
##'             control = list(pch = 20, col = "blue", PF.pch = 17,
##'                            PF.points.col = "blue", PF.line.col = "blue"))
##' text(omEGO1$values[,1], omEGO1$values[,2], labels = 1:nsteps, pos = 3, col = "blue")
##' plotGPareto(list = omEGO2, add = TRUE,
##'             control = list(pch = 20, col = "green", PF.pch = 17,
##'                            PF.points.col = "green", PF.line.col = "green"))
##' text(omEGO2$values[,1], omEGO2$values[,2], labels = 1:nsteps, pos = 3, col = "green")
##' plotGPareto(list = omEGO3, add = TRUE,
##'             control = list(pch = 20, col = "red", PF.pch = 17,
##'                            PF.points.col = "red", PF.line.col = "red"))
##' text(omEGO3$values[,1], omEGO3$values[,2], labels = 1:nsteps, pos = 3, col = "red") 
##' plotGPareto(list = omEGO4, add = TRUE,
##'             control = list(pch = 20, col = "orange", PF.pch = 17,
##'                            PF.points.col = "orange", PF.line.col = "orange"))
##' text(omEGO4$values[,1], omEGO4$values[,2], labels = 1:nsteps, pos = 3, col = "orange")
##' points(response.grid[,1], response.grid[,2], col = "black", pch = 20)
##' legend("topright", c("EHI", "SMS", "SUR", "EMI"), col = c("blue", "green", "red", "orange"),
##'  pch = rep(17,4))
##'  
##'  
##' # Post-processing
##' plotGPareto(omEGO1, UQ_PF = T, UQ_PS = T, UQ_dens = T)
##'   
##' }

GParetoptim <- function (model, fn, cheapfn=NULL, crit="SMS", nsteps, lower, upper, type="UK", cov.reestim=TRUE,
                         critcontrol = NULL,
                         optimcontrol = list(method="genoud", threshold = 1e-5, distance = "euclidean", notrace = FALSE), ...){
  ##########################################################################################
  # Inputs :
  # model: list of 2 or 3 models
  # fn: objective function, returns 2 or 3 objectives
  # nsteps: number of iterations
  # lower, upper: design region
  # optimcontrol, type, CovReEstimate: parameters as in DiceOptim & KrigInv
  ##########################################################################################
  n     <- nrow(model[[1]]@X)
  d     <- model[[1]]@d
  
  fn <- match.fun(fn)
  
  if(is.null(optimcontrol$method))
    optimcontrol$method <- "genoud"
  
  # Build fastfun if necessary
  if (!is.null(cheapfn)) {
    cheapfn <- match.fun(cheapfn)
    fastobs <- apply(model[[1]]@X, 1, cheapfn)
    fastmod <- fastfun(fn = cheapfn, design = model[[1]]@X, response = fastobs)
    model[[length(model)+1]] <- fastmod
  }else{
    Y.new.cheap = NULL
  }
  n.obj <- length(model)
  
  if (length(model) < 2 ){
    cat("Error in model definition: 'model' must be a list of 2 or more km models \n")
    return(NULL)
  }
  
  if (length(model) >= 3 && crit=="EHI"){
    cat("Analytical Hypervolume EI only works with 2 objectives; SAA approximation used. \n")
  }
  if(length(model) > 3 && crit == "SUR"){
    cat("SUR is available for 2 or 3 objectives \n")
  }
  
  # Regroup all observations
  observations <- c()
  for (i in 1:n.obj) observations <- cbind(observations, model[[i]]@y)
  
  if(is.null(optimcontrol$notrace)){
    notrace <- FALSE
  }else{
    notrace <- optimcontrol$notrace
  }
  
  if(!notrace){
    cat("----------------------------\n")
    cat("Starting optimization with : \n The criterion", crit, "\n The solver",  optimcontrol$method, "\n")
    cat("----------------------------\n")
    cat("Ite / Crit / New x / New y \n")
  }
  paretoFront <- t(nondominated_points(t(observations)))
  noRef <- FALSE #TRUE if a refPoint is provided
  if(is.null(critcontrol$refPoint) & (crit == "SMS" | crit =="EHI")){
    noRef <- TRUE
    estimatedRef <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj) ### May be changed
    critcontrol$refPoint <- estimatedRef 
    if((crit == "SMS" | crit =="EHI") & !notrace)
      cat("No refPoint provided, ", signif(critcontrol$refPoint, 3), "used \n")
  }
  
  ## Change of seed for optimization unless set for genoud (for now by setting default values)
  ## NOTE: might be interesting to modify the defaults in crit_optimizer
  if(is.null(optimcontrol$unif.seed))
    changeSeed1 <- TRUE
  if(is.null(optimcontrol$int.seed))
    changeSeed2 <- TRUE
  
  #### Main loop starts here ################
  for (i in 1:nsteps) {
    
    ## Compute current Pareto front
    paretoFront <- t(nondominated_points(t(observations)))
    n.pareto    <- nrow(paretoFront)
    
    ## Change the refPoint if none is provided, if necessary
    if(noRef & (crit == "SMS" | crit =="EHI") & !notrace){
      if(any(estimatedRef != matrix(apply(paretoFront, 2, max) + 1, 1, n.obj))){
        estimatedRef <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj) ### May be changed
        critcontrol$refPoint <- estimatedRef
        cat("refPoint changed, ", signif(critcontrol$refPoint, 3), "used \n")
      }
    }
      
    
    ## Change the seeds for genoud to avoid selecting always the same initial values
    if(optimcontrol$method == "genoud" & changeSeed1){
      optimcontrol$unif.seed <- sample.int(1e9, 1)
    }
    
    if(optimcontrol$method == "genoud" & changeSeed2){
      optimcontrol$int.seed <- sample.int(1e9, 1)
    }
    

    # observations removed (could be reintegrated)
    sol <- try(crit_optimizer(crit=crit, model=model, lower=lower, upper=upper, 
                              optimcontrol=optimcontrol, type=type, paretoFront=paretoFront, 
                              critcontrol=critcontrol, nsteps.remaining=nsteps-i))
    
    if (typeof(sol) == "character") {
      if(!notrace){
        cat("Unable to maximize criterion at iteration ", i, "- optimization stopped \n")
        cat("Last model returned \n")
      }
      
      par <- values <- c()
      if (i > 1) {
        par <- model[[1]]@X[(n+1):model[[1]]@n,, drop=FALSE]
        # par <- model[[1]]@X[(n+1):model[[1]]@n,, drop=FALSE]
      }
      return(list(par=par, values=values, nsteps = i-1, lastmodel = model))
    }
    
    ## Check if optimization do not return already known point
    if(checkPredict(sol$par, model, type = type, distance = critcontrol$distance, threshold = critcontrol$threshold)){
      if(!notrace)
        cat("Optimization failed, so a random point is selected (consider increasing the inner optimization budget).\n")
      sol$par <- matrix(runif(d), nrow = 1)
    }
      
    
    ## Update
    X.new <- matrix(as.numeric(sol$par), nrow=1, ncol=d)
    Y.new <- try(fn(as.numeric(sol$par), ...))
    if (!is.null(cheapfn)) {
      Y.new.cheap <- try(cheapfn(as.numeric(sol$par)))
    }
    
    if (typeof(Y.new) == "character" || (!is.null(cheapfn) && typeof(Y.new.cheap) == "character")) {
      if(!notrace){
        cat("Unable to compute objective function at iteration ", i, "- optimization stopped \n")
        cat("Problem occured for the design: ", X.new, "\n")
        cat("Last model returned \n")
      }
      
      par <- values <- c()
      if (i > 1) {
        # par <- model[[1]]@X[(n+1):model[[1]]@n,, drop=FALSE]
        par <- model[[1]]@X[(n+1):model[[1]]@n,, drop=FALSE]
      }
      return(list(par=par, values=values, nsteps = i-1, lastmodel = model))
    }
    Y.new <- c(Y.new, Y.new.cheap)
    if(!notrace) cat( i, "/", signif(sol$val,3), "/", signif(X.new,3), "/", signif(Y.new,3), "\n", sep = "\t")
    
    # Remove new observation from integration points if discrete case is used
    if (optimcontrol$method=="discrete") {
      optimcontrol$candidate.points <- optimcontrol$candidate.points[-sol$index,,drop=FALSE]
      if (crit=="SUR") { 
        critcontrol$integration.points <- critcontrol$integration.points[-sol$index,,drop=FALSE]
      }
    }
    #     print(model[[1]]@covariance)
    # Update models
    observations <- rbind(observations, Y.new)
    newmodel <- model
    for (j in 1:n.obj) {
      newmodel[[j]] <- try(update(object = model[[j]], newX = X.new, newy=Y.new[j], newX.alreadyExist=FALSE,
                                  cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE))), silent = TRUE)
      if (typeof(newmodel[[j]]) == "character" && cov.reestim) {
        cat("Error in hyperparameter estimation - old hyperparameter values used instead for model ", j, "\n")
        newmodel[[j]] <- try(update(object = model[[j]], newX = X.new, newy=Y.new[j], newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
      }
      if (typeof(newmodel[[j]]) == "character") {
        cat("Unable to udpate kriging model ", j, " at iteration", i-1, "- optimization stopped \n")
        cat("lastmodel ", j, " is the model at iteration", i-1, "\n")
        cat("par and values contain the ",i , "th observation \n \n")
        if (i > 1) allX.new <- rbind(model[[1]]@X[(n+1):(n+i-1),, drop=FALSE], X.new)
        return(list(
          par    = allX.new,
          values = observations[(n+1):(n+i),, drop=FALSE],
          nsteps = i, 
          lastmodel = model))
      } else {
        model[[j]] <- newmodel[[j]]
      }
    }
    
  }
  
  if(!notrace) cat("\n")
  #### End of main loop ################
  return(list(
    par=model[[1]]@X[(n+1):(n+nsteps),, drop=FALSE], 
    values=observations[(n+1):(n+nsteps),, drop=FALSE], 
    nsteps=nsteps, 
    lastmodel=model))
}
