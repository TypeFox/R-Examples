##' User-friendly wrapper of the function \code{\link[GPareto]{GParetoptim}}.
##' Generates initial DOEs and kriging models (objects of class \code{\link[DiceKriging]{km}}), 
##' and executes \code{nsteps} iterations of multiobjective EGO methods.
##' @title EGO algorithm for multiobjective optimization
##' @details Does not require specific knowledge on kriging models (objects of class \code{\link[DiceKriging]{km}}).\cr
##' 
##' The problem considered is of the form: \eqn{min f(x) = f_1(x), ..., f_p(x)}. 
##' The \code{control} argument is a list that can supply any of the following optional components: \cr
##' \itemize{
##' \item \code{method}{: choice of multiobjective improvement function: "\code{SMS}", "\code{EHI}", "\code{EMI}" or "\code{SUR}" 
##' (see \code{\link[GPareto]{crit_SMS}}, \code{\link[GPareto]{crit_EHI}}, \code{\link[GPareto]{crit_EMI}}, \code{\link[GPareto]{crit_SUR}}),}
##' \item \code{trace}{:  if positive, tracing information on the progress of the optimization is produced,}
##' \item \code{inneroptim}{: choice of the inner optimization algorithm: "\code{genoud}", "\code{pso}" or "\code{random}"
##'  (see \code{\link[rgenoud]{genoud}} and \code{\link[pso]{psoptim}}),}
##' \item \code{maxit}{: maximum number of iterations of the inner loop,}
##' \item \code{seed}{: to fix the random variable generator,}
##' \item \code{refPoint}{: reference point for hypervolume computations (for "\code{SMS}" and "\code{EHI}" methods).}
##' }
##' For additional details or other possible arguments, see \code{\link[GPareto]{GParetoptim}}.\cr
##' 
##' Display of results and various post-processings are available with \code{\link[GPareto]{plotGPareto}}.  
##' 
##' @param fn the multi-objective function to be minimized (vectorial output), found by a call to \code{\link[base]{match.fun}},
##' @param cheapfn optional additional fast-to-evaluate objective function (handled next with class \code{\link[GPareto]{fastfun}}), 
##' which does not need a kriging model, handled by a call to \code{\link[base]{match.fun}},
##' @param budget total number of calls to the objective function,
##' @param lower vector of lower bounds for the variables to be optimized over,
##' @param upper vector of upper bounds for the variables to be optimized over,
##' @param par initial design of experiments. If not provided, \code{par} is taken as a maximin LHD with budget/3 points,
##' @param value initial set of objective observations \eqn{f(par)}. Computed if not provided.
##' Not that \code{value} may NOT contain any \code{cheapfn} value,
##' @param control an optional list of control parameters. See "Details",
##' @param ... additional parameters to be given to the objective \code{fn}.
##' @export
##' @return
##' A list with components:
##' \itemize{
##' \item{\code{par}}{: all the non-dominated points found,}
##' \item{\code{value}}{: the matrix of objective values at the points given in \code{par},}
##' \item{\code{history}}{: a list containing all the points visited by the algorithm (\code{X}) and their corresponding objectives (\code{y}),}
##' \item{\code{model}}{: a list of objects of class \code{\link[DiceKriging]{km}}, corresponding to the last kriging models fitted.}
##' }
##' 
##' 
##' @author
##' Victor Picheny (INRA, Castanet-Tolosan, France)
##' 
##' Mickael Binois (Mines Saint-Etienne/Renault, France)
##' 
##' @references 
##' M. T. Emmerich, A. H. Deutz, J. W. Klinkenberg (2011), Hypervolume-based expected improvement: Monotonicity properties and exact computation,
##' \emph{Evolutionary Computation (CEC)}, 2147-2154. \cr \cr
##' V. Picheny (2015), Multiobjective optimization using Gaussian process emulators via stepwise uncertainty reduction, 
##' \emph{Statistics and Computing}, 25(6), 1265-1280. \cr \cr
##' T. Wagner, M. Emmerich, A. Deutz, W. Ponweiser (2010), On expected-improvement criteria for model-based multi-objective optimization.   
##' \emph{Parallel Problem Solving from Nature}, 718-727, Springer, Berlin. \cr \cr
##' J. D. Svenson (2011), \emph{Computer Experiments: Multiobjective Optimization and Sensitivity Analysis}, Ohio State university, PhD thesis. 
##' 
##' @importFrom DiceDesign maximinESE_LHS lhsDesign
##' 
##' @examples
##' #---------------------------------------------------------------------------
##' # 2D objective function, 3 cases
##' #---------------------------------------------------------------------------
##' \dontrun{
##' set.seed(25468)
##' n_var <- 2 
##' fname <- ZDT3
##' lower <- rep(0, n_var)
##' upper <- rep(1, n_var)
##' 
##' #---------------------------------------------------------------------------
##' # 1- Expected Hypervolume Improvement optimization, using pso
##' #---------------------------------------------------------------------------
##' res <- easyGParetoptim(fn=fname, lower=lower, upper=upper, budget=15, 
##'                    control=list(method="EHI", inneroptim="pso", maxit=20))
##' par(mfrow=c(1,2))
##' plotGPareto(res)
##' title("Pareto Front")
##' plot(res$history$X, main="Pareto set", col = "red", pch = 20)
##' points(res$par, col="blue", pch = 17)
##' 
##' #---------------------------------------------------------------------------
##' # 2- SMS Improvement optimization using random search, with initial DOE given
##' #---------------------------------------------------------------------------
##' library(DiceDesign)
##' design.init   <- maximinESE_LHS(lhsDesign(10, n_var, seed = 42)$design)$design
##' response.init <- t(apply(design.init, 1, fname))
##' res <- easyGParetoptim(fn=fname, par=design.init, value=response.init, lower=lower, upper=upper, 
##'                        budget=15, control=list(method="SMS", inneroptim="random", maxit=100))
##' par(mfrow=c(1,2))
##' plotGPareto(res)
##' title("Pareto Front")
##' plot(res$history$X, main="Pareto set", col = "red", pch = 20)
##' points(res$par, col="blue", pch = 17)
##' 
##' #---------------------------------------------------------------------------
##' # 3- Stepwise Uncertainty Reduction optimization, with one fast objective function
##' #---------------------------------------------------------------------------
##' fname <- camelback
##' cheapfn <- function(x) {
##' if (is.null(dim(x))) return(-sum(x))
##' else return(-rowSums(x))
##' }
##' res <- easyGParetoptim(fn=fname, cheapfn=cheapfn, lower=lower, upper=upper, budget=15, 
##'                    control=list(method="SUR", inneroptim="pso", maxit=20))
##' par(mfrow=c(1,2))
##' plotGPareto(res)
##' title("Pareto Front")
##' plot(res$history$X, main="Pareto set", col = "red", pch = 20)
##' points(res$par, col="blue", pch = 17)
##' }

easyGParetoptim <- function (fn, cheapfn = NULL, budget, lower, upper, par=NULL, value=NULL,  
                             control=list(method="SMS", trace=1, inneroptim="pso", maxit=100, seed=42), ...) {
  
  if (is.null(control$method)) control$method <- "SMS"
  if (is.null(control$trace)) control$trace   <- 1
  if (is.null(control$inneroptim))  control$inneroptim <- "pso"
  if (is.null(control$maxit)) control$maxit   <- 100
  if (is.null(control$seed)) control$seed   <- 42
  if (is.null(control$refPoint)) control$refPoint   <- NULL
  if (is.null(control$distance)) control$distance <- "euclidean"
  
  fn <- match.fun(fn)
  d <- length(lower)
  
  if (length(lower) != length(upper)) {
    cat("Bound values lower and upper are not consistent. Both should be vectors of size d.")
    return(0)
  }
  
  if (!is.null(par)) {
    design.init <- data.frame(x=par) #, ncol=d, byrow=TRUE)
    temp <- dim(design.init)
    
    if (temp[2] != d) {
      cat("Bound values (lower and upper) and initial DOE (par) are not consistent. \n 
          lower and upper should be vectors of size d and \n 
          par either a matrix with d columns or a data frame of d variables.")
      return(0)
    }
    n.init <- temp[1]
  } else {
    n.init <- max(4*d, round(budget/3))
    design.init <- data.frame(x=lower + (upper-lower)*maximinESE_LHS(lhsDesign(n.init, d, seed=control$seed)$design)$design)
  }
  #----------------------------------------
  if (!is.null(par) && !is.null(value)) {
    obs.init <- as.matrix(value)
    if (nrow(obs.init) != n.init) {
      cat("Initial DOE (par) and objective (value) are not consistent.")
      return(0)
    }
  } else {
    obs.init <- as.matrix(apply(design.init, 1, fn, ...))
    if (nrow(obs.init)==1) obs.init <- t(obs.init)
    if (nrow(obs.init)!=n.init) obs.init <- t(obs.init)
  }
  n.obj <- ncol(obs.init)
  #----------------------------------------
  if (!is.null(par) && !is.null(value)) {
    n.ite <- budget
  } else {
    n.ite <- budget - n.init
  }
  #----------------------------------------
  model <- vector("list", n.obj)
  for (j in 1:(n.obj)) {
    model[[j]] <- km(~., design = design.init, response = obs.init[,j], control=list(trace=FALSE), lower=rep(.1,d), upper=rep(1,d))
  }
  #----------------------------------------
  critcontrol <- NULL

  if (control$inneroptim=="genoud") optimcontrol <- list(method="genoud", max.generations=control$maxit, notrace = !control$trace>0)
  if (control$inneroptim=="pso")    optimcontrol <- list(method="pso", maxit=control$maxit, notrace = !control$trace>0)
  if (control$inneroptim=="random") {
    candidate.points <- matrix(rep(lower,control$maxit)+rep(upper-lower,control$maxit)*runif(control$maxit*d),byrow=T,ncol=d)
    optimcontrol = list(method="discrete", candidate.points=candidate.points, notrace = !control$trace>0)
  }
  
  if (control$method=="SUR") {
    critcontrol$distrib <- "SUR"
    critcontrol$n.points <- 50*d
  } else if (control$method %in% c("SMS", "EHI")) {
    critcontrol$refPoint <- control$refPoint
  }
  
  if(!is.null(control$distance)) critcontrol$distance <- control$distance
  if(!is.null(control$threshold)) critcontrol$threshold <- control$threshold

  omEGO <- GParetoptim(model = model, fn = fn, cheapfn = cheapfn,
                       crit = control$method, nsteps=n.ite, lower=lower, upper=upper, cov.reestim=TRUE, 
                       optimcontrol=optimcontrol, critcontrol = critcontrol, ...)
  
  allX <- omEGO$lastmodel[[1]]@X
  
  ally <- c()
  ## the last row is given by omEGOvalues, since some models may not have been updated  
  for (i in 1:length(omEGO$lastmodel)) ally <- cbind(ally,
                                                     omEGO$lastmodel[[i]]@y[1:(omEGO$lastmodel[[1]]@n-1)])
  ally <- rbind(ally, omEGO$values[nrow(omEGO$values),])

  # Compute current best
  value <- t(nondominated_points(t(ally)))
  par   <- allX[!is_dominated(t(ally)),]
  
  return(list(par=par, value = value, history=list(X=allX, y=ally, model=omEGO$lastmodel)))
  }