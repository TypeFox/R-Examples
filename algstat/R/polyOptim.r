#' Polynomial Optimization
#'
#' Find the collection of critical points of a multivariate polynomial unconstrained or constrained to an affine variety (algebraic set; solution set of multivariate polynomials).
#' 
#' @param objective the objective polynomial (as a character or mpoly)
#' @param constraints (as a character or mpoly/mpolyList)
#' @param varOrder variable order (see examples)
#' @param ... stuff to pass to bertini
#' @return an object of class bertini
#' @export polyOptim
#' @examples
#' \dontrun{
#' 
#' # unconstrained optimization of polynomial functions is available
#' polyOptim("x^2") 
#' polyOptim("-x^2") 
#' polyOptim("-(x - 2)^2") 
#' polyOptim("-(x^2 + y^2)") 
#' polyOptim("-(x^2 + (y - 2)^2)") 
#'
#' polyOptim("(x - 1) (x - 2) (x - 3)") # fix global labeling
#'
#'
#' # constrained optimization over the affine varieties is also available
#' # (affine variety = solution set of polynomial equations)
#'
#' # find the critical points of the plane f(x,y) = x + y
#' # over the unit circle x^2 + y^2 = 1
#' polyOptim("x + y", "x^2 + y^2 = 1")
#'
#' # you can specify them as a combo of mpoly, mpolyList, and characters
#' o <- mp("x + y")
#' c <- "x^2 + y^2 = 1"
#' polyOptim(o, c)
#'
#' c <- mp("x^2 + y^2 - 1")
#' polyOptim(o, c)
#'
#' out <- polyOptim("x + y", c)
#' str(out)
#' 
#' # another example, note the solutions are computed over the complex numbers
#' polyOptim("x^2 y", "x^2 + y^2 = 3")
#' # solutions: (+-sqrt(2), +-1) and (0, +-sqrt(3))
#'
#'
#' 
#'
#' }
#'
polyOptim <- function(objective, constraints, varOrder, ...){

  optimizationType <- "unconstrained"
  
  if(!missing(constraints)){
  	
  	optimizationType <- "constrained"
    
    if(is.character(constraints)){
      if(all(str_detect(constraints, "=="))){
        split <- strsplit(constraints, " == ")
        lhs <- sapply(split, function(x) x[1])
        rhs <- sapply(split, function(x) x[2])
        constraints <- mp(lhs) - mp(rhs)
      } else if(all(str_detect(constraints, "="))){
        split <- strsplit(constraints, " = ")      
        lhs <- sapply(split, function(x) x[1])
        rhs <- sapply(split, function(x) x[2])                
        constraints <- mp(lhs) - mp(rhs)      
      } else {
        constraints <- mp(constraints)
      }
    }
        
    if(is.mpoly(constraints)){
      constraints <- list(constraints)
      class(constraints) <- "mpolyList"
    }
  
    # add lagrange multipliers
    lams <- paste0("l", length(constraints))
    mults <- mp( lams )
    if(is.mpoly(mults)){
      mults <- list(mults)
      class(mults) <- "mpolyList"
    }
    lagrangeConstraints <- mults * constraints
  } else {
    lams <- NULL
  }


  if(is.character(objective)) objective <- mp(objective)
  stopifnot(is.mpoly(objective))  


  # sort out variables
  objectiveVars <- vars(objective)
  nVars <- length(objectiveVars)
  nLagrangeMults <- length(lams)
  vars <- c(objectiveVars, lams)
  
  if(!missing(varOrder) && 
    !all(sort(objectiveVars) == sort(varOrder))
  ){
    stop(
      "if varOrder is provided, it must contain all of the variables.", 
      call.=F
    )
  }
  if(!missing(varOrder)) vars <- varOrder

  if(optimizationType == "unconstrained"){
  	grad <- deriv(objective, var = vars)
  }
  
  if(optimizationType == "constrained"){
    lagrangian <- objective + Reduce("+", lagrangeConstraints)
  	grad <- deriv(lagrangian, var = vars)
  }  

  out <- polySolve(grad, varOrder = vars)

  # add optim related stuff to the output
  out$variables <- list(vars = objectiveVars, lams = lams)
  f <- suppressMessages(as.function(objective, varorder = objectiveVars))
  real_optima <- as.data.frame(out$real_finite_solutions)
  real_optima$value <- apply(real_optima[,1:nVars,drop=FALSE], 1, f)
  real_optima <- real_optima[order(real_optima$value, decreasing = TRUE),]
  real_optima$optima <- ""
  real_optima$optima[which.max(real_optima$value)] <- "global max"
  real_optima$optima[which.min(real_optima$value)] <- "global min"  
  if(optimizationType == "unconstrained") real_optima$optima <- ""  
  out$real_optima <- real_optima

  class(out) <- c("polyOptim", "bertini")

  out
}