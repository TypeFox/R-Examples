##' Modification of an R function to be used as with methods \code{predict} and \code{update} (similar to a \code{\link[DiceKriging]{km}} object). 
##' It creates an S4 object which contains the values corresponding to evaluations of other costly observations.
##' It is useful when an objective can be evaluated fast.
##' @title Fastfun function
##' @param fn the evaluator function, found by a call to \code{\link[base]{match.fun}},
##' @param design a data frame representing the design of experiments.
##' The ith row contains the values of the d input variables corresponding to the ith evaluation.
##' @param response optional vector (or 1-column matrix or data frame) containing the values of the 1-dimensional output given by the objective function at the design points.
##' @import methods 
##' @export
##' @examples
##' ########################################################
##' ## Example with a fast to evaluate objective
##' ########################################################
##' \dontrun{
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' d <- 2 
##' 
##' fname <- P1
##' n.grid <- 21
##' nappr <- 11 
##' design.grid <- maximinESE_LHS(lhsDesign(nappr, d, seed = 42)$design)$design
##' response.grid <- t(apply(design.grid, 1, fname))
##' Front_Pareto <- t(nondominated_points(t(response.grid)))
##' 
##' mf1 <- km(~., design = design.grid, response = response.grid[,1])
##' mf2 <- km(~., design = design.grid, response = response.grid[,2])
##' model <- list(mf1, mf2)
##' 
##' nsteps <- 5 
##' lower <- rep(0, d)
##' upper <- rep(1, d)
##' 
##' # Optimization reference: SMS with discrete search
##' optimcontrol <- list(method = "pso")
##' omEGO1 <- GParetoptim(model = model, fn = fname, crit = "SMS", nsteps = nsteps,
##'                      lower = lower, upper = upper, optimcontrol = optimcontrol)
##' print(omEGO1$par)
##' print(omEGO1$values)
##' plot(response.grid, xlim = c(0,300), ylim = c(-40,0), pch = 17, col = "blue")
##' points(omEGO1$values, pch = 20, col ="green") 
##' 
##' # Optimization with fastfun: SMS with discrete search
##' # Separation of the problem P1 in two objectives: 
##' # the first one to be kriged, the second one with fastobj
##' f1 <-   function(x){
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1) 
##'   b1 <- 15*x[,1] - 5
##'   b2 <- 15*x[,2]
##'   return(  (b2 - 5.1*(b1/(2*pi))^2 + 5/pi*b1 - 6)^2 +10*((1 - 1/(8*pi))*cos(b1) + 1))
##' }
##' 
##' f2 <-   function(x){
##'   if(is.null(dim(x))) x <- matrix(x, nrow = 1) 
##'   b1<-15*x[,1] - 5
##'   b2<-15*x[,2]
##'   return(-sqrt((10.5 - b1)*(b1 + 5.5)*(b2 + 0.5))
##'          - 1/30*(b2 - 5.1*(b1/(2*pi))^2 - 6)^2
##'          - 1/3*((1 - 1/(8*pi))*cos(b1) + 1)) 
##' }
##' 
##' optimcontrol <- list(method = "pso")
##' model2 <- list(mf1)
##' omEGO2 <- GParetoptim(model = model2, fn = f1, cheapfn = f2, crit = "SMS", nsteps = nsteps,
##'                      lower = lower, upper = upper, optimcontrol = optimcontrol)
##' print(omEGO2$par)
##' print(omEGO2$values)
##' 
##' points(omEGO2$values, col = "red", pch = 15)
##' }
##' @return An object of class  \code{\link[GPareto]{fastfun-class}}.
`fastfun` <-
  function(fn, design, response = NULL) {
    
    model <- new("fastfun")
    
    data <- data.frame(design)
    
    X <- as.matrix(data)
    
    fn <- match.fun(fn)
    
    if(is.null(response)){
      response <- apply(design, 1, fn)
    }
    y <- matrix(response, ncol=1)
    model@X <- X
    model@y <- y
    model@d <- ncol(X)
    model@n <- nrow(X)
    model@fun <- fn
    
    return(model)
  }

## ----------------
## CLASS definition
## ----------------

## fastfun Class
##' Class for fast to compute objective.
##' @slot d spatial dimension,
##' @slot n observations number,
##' @slot X the design of experiments, size \code{n x d},
##' @slot y  the observations, size \code{n x 1},
##' @slot fun the evaluator function.
##' @section Objects from the Class : To create a \code{fastfun} object, use \code{\link[GPareto]{fastfun}}. See also this function for more details and examples.
##' @export
setClass("fastfun", 		
         representation( 
           d = "integer",          ## spatial dimension
           n = "integer",          ## observations number
           ## data
           X = "matrix",           ## the design of experiments, size nxd
           y = "matrix",           ## the observations, size nx1
           fun = "function"        ## the evaluator
           
         )
)

# if(!isGeneric("show")) {
#   setGeneric(name = "show",
#              def = function(object) standardGeneric("show")
#   )
# }
# 
# setMethod("show", "fastfun", 
#           function(object){
#             show.fastfun(object)		
#           }
# )


##*****************************************************************************
##                        P R E D I C T  METHOD
##*****************************************************************************


predict.fastfun <- function(object, newdata, cov.compute = TRUE, light.return = FALSE, ...) {
  res <- list(trend = NULL, mean = apply(newdata, 1, object@fun), sd = rep(0, nrow(newdata)))
  if(!light.return & cov.compute){
    res$cov <- matrix(0, nrow(newdata), nrow(newdata))
  }
  res$trend <- res$mean
  res$lower95 <- res$mean
  res$upper95 <- res$mean
  
  return(res)
  
}



if(!isGeneric("predict")) {
  setGeneric(name = "predict",
             def = function(object, ...) standardGeneric("predict")
  )
}

##' @describeIn fastfun Predict(by evaluating \code{fun}) the result at a new observation.
##' @param newdata Matrix of the new location for the design
setMethod("predict", "fastfun", 
          function(object, newdata, ...) {
            predict.fastfun(object = object, newdata = newdata, ...)
          }
)


##****************************************************************************
##			     update  METHOD
##****************************************************************************

update.fastfun <- function(object, newX, newy, ...){
    object@X <- rbind(object@X, matrix(newX, ncol=ncol(object@X), byrow=TRUE))
    object@y <- rbind(object@y, matrix(newy, ncol=1, byrow=TRUE))
    
    object@n <- as.integer(object@n + 1)
    
    return(object)
}

if(!isGeneric("update")) {
  setGeneric(name = "update",
             def = function(object, ...) standardGeneric("update")
  )
}

##' @param object \code{\link[GPareto]{fastfun}} object
##' @param newX Matrix of the new location for the design
##' @param newy Matrix of the responses at \code{newX}
##' @param ... further arguments (not used)
##' @describeIn fastfun Update the \code{X} and \code{y} slots with a new design and observation.
##' @keywords internal
setMethod("update", "fastfun", 
          function(object, newX, newy, ...) {
            update.fastfun(object = object, newX = newX, newy = newy, ...) 
          }
)


##****************************************************************************
##			     simulate  METHOD
##****************************************************************************

simulate.fastfun <- function(object, nsim=1, seed=NULL, newdata=NULL, 
                             cond=FALSE, nugget.sim=0, checkNames=TRUE, ...){
  if(is.null(newdata)){
    simulation <- as.vector(object@y)
  }else{
    simulation <- as.vector(apply(newdata, 1, object@fun))
  }
  return(matrix(simulation, nrow = nsim, ncol = length(simulation), byrow = T))
}

if(!isGeneric("simulate")) {
  setGeneric(name = "simulate",
             def = function(object, ...) standardGeneric("simulate")
  )
}

##' @param object \code{\link[GPareto]{fastfun}} object
##' @param nsim an optional number specifying the number of response vectors to simulate. Default is 1.
##' @param seed usual seed argument of method simulate. Not used.
##' @param newdata an optional vector, matrix or data frame containing the points where to perform predictions.
##'  Default is \code{NULL}: simulation is performed at design points specified in \code{object}.
##' @param cond an optional boolean indicating the type of simulations. Not used.
##' @param nugget.sim	an optional number corresponding to a numerical nugget effect. Not used.
##' @param checkNames an optional boolean. Not used.
##' @param ... further arguments (not used)
##' @describeIn fastfun Simulate responses values (for compatibility with methods using \code{\link[DiceKriging]{simulate}})
##' @keywords internal
setMethod("simulate", "fastfun", 
          function(object, nsim, seed, newdata, cond, nugget.sim, checkNames, ...) {
            simulate.fastfun(object = object, nsim = nsim, seed = seed, newdata = newdata, 
                             cond = cond, nugget.sim = nugget.sim, checkNames = checkNames, ...) 
          }
)
