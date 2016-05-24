library(slam)
library(compiler)

#' Create a new LP or MIP
#' 
#' \code{newProblem} creates a new and empty linear program (LP)
#' or mixed integer program (MIP).
#' 
#' @param max \code{TRUE} (default) for a maximization problem, \code{FALSE} for a minization problem
#' @return An object of class \code{lpmodeler}.
#' @author Cyrille Szymanski <cnszym at gmail.com>
#' @seealso TODO
#' @export
#' @import slam
#' @examples
#' p <- newProblem()
newProblem <- function(max=T)
{
  prob <- list()
  prob$max <- max
  prob$obj <- c()
  prob$mat <- simple_triplet_zero_matrix(nrow=0,ncol=0)
  prob$dir <- c()
  prob$rhs <- c()
  prob$types <- c()
  prob$vars <- new.env()
  prob$ctrs <- new.env()
  class(prob) <- "lpmodeler"
  prob
}

#' Print a LP or MIP problem
#' 
#' Prints general information about a linear program (LP) or 
#' mixed integer program (MIP) represented by an object of class
#' \code{lpmodeler}.
#' 
#' @param x an object of class \code{lpmodeler}
#' @param ... further arguments passed to or from other methods
#' @author Cyrille Szymanski <cnszym at gmail.com>
#' @seealso TODO
#' @export
#' @examples
#' p <- newProblem()
#' print(p)
print.lpmodeler <- function(x, ...) {
  cat("lpmodeler: ",
      NCOL(x$mat), " variables and ", 
      NROW(x$mat), " constraints.\n",
      sep="")
  invisible(x)
}

#' Check the consistency of the dimensions of a LP or MIP
#' 
#' \code{checkDims} checks the consistency of the dimensions of
#' the matrices and vectors of a linear program (LP) or a mixed
#' integer program (MIP) reprensented by an object of class
#' \code{lpmodeler}.
#' 
#' @param p an object of class \code{lpmodeler}
#' @author Cyrille Szymanski <cnszym at gmail.com>
#' @seealso TODO
#' @export
checkDims <- function(p)
{
  if( length(p$obj)!=NCOL(p$mat) ) stop("incorrect dimensions: obj")
  if( length(p$types)!=NCOL(p$mat) ) stop("incorrect dimensions: types")
  if( length(p$dir)!=NROW(p$mat) ) stop("incorrect dimensions: dir")
  if( length(p$rhs)!=NROW(p$mat) ) stop("incorrect dimensions: rhs")
}

#' Add a new variable to a LP or MIP
#' 
#' \code{addVariable} creates a new variable (continuous, integer
#' or binary) and adds it to a linear program (LP) or mixed integer
#' program (MIP) represented by an object of class \code{lpmodeler}.
#' 
#' @details
#' TODO
#' 
#' @param p an object of class \code{lpmodeler}
#' @param t type of the variable to create, C = continuous (default), I = integer, B = binary
#' @param o numeric value representing the coefficient of the variable in the objective function (objective point), 0 by default
#' @param name an optional string to name the new variable
#' @return An object of class \code{lpmodeler}.
#' @author Cyrille Szymanski <cnszym at gmail.com>
#' @seealso TODO
#' @export
#' @examples
#' p <- newProblem()
#' 
#' # add an integer variable called x to the
#' # problem, and set its coefficient in the
#' # objective function to a value of 5.
#' p <- addVariable(p, "I", 5, "x")
addVariable <- cmpfun(function(p,t=c("C","I","B"),o=0,name=NULL)
{
  t <- t[1]
  if( !t %in% c("C","I","B")) stop("invalid type, should be equal to C, I or B")
  p$obj <- c(p$obj,o)
  p$mat <- cbind(p$mat,simple_triplet_zero_matrix(nrow=NROW(p$mat),ncol=1))
  p$types <- c(p$types,t)
  if( !is.null(name) ) 
    assign(name, NCOL(p$mat), envir=p$vars) # TODO vérifier qu'elle n'existe pas déjà
  p
})

#' Add a new constraint to a LP or MIP
#' 
#' \code{addConstraint} creates a new constraint (<, >, <=, >=, ==)
#' and adds it to a linear program (LP) or mixed integer
#' program (MIP) represented by an object of class \code{lpmodeler}.
#' 
#' @details
#' TODO
#' 
#' @param p an object of class \code{lpmodeler}
#' @param sense sense of the constraint (\code{<}, \code{>}, \code{<=}, \code{>=}, \code{==} or \code{!=})
#' @param rhs right hand side of the constraint
#' @param coefs optional coefficients of the variables in the left hand side of the constraint
#' @param name an optional string to name the new constraint
#' @return An object of class \code{lpmodeler}.
#' @author Cyrille Szymanski <cnszym at gmail.com>
#' @seealso TODO
#' @export
#' @examples
#' p <- newProblem()
#' 
#' # add variables x and y
#' p <- addVariable(p, "C", 5, "x")
#' p <- addVariable(p, "C", 4, "y")
#' 
#' # add the constraint: x + 2y >= 5
#' p <- addConstraint(p, ">=", 5, c(1, 2), name = "x+2y greater or equal than 5")
#' 
#' # add the empty constraint: <= 10
#' p <- addConstraint(p, "<=", 10, name = "less or equal than 10")
addConstraint <- cmpfun(function(p,sense,rhs,coefs=NULL,name=NULL)
{
  if( !sense %in% c("<=",">=","<",">","==","!=") ) stop("invalid sense,should be equal to <, >, <=, >= or ==")
  if( !is.null(coefs) ) {
    if( length(coefs) != NCOL(p$mat) ) stop("invalid dimensions, problem has ",NCOL(p$mat)," variables but ",length(coefs)," coefficients supplied")
    p$mat <- rbind(p$mat,as.simple_triplet_matrix(matrix(coefs,nrow=1,ncol=NCOL(p$mat))))
  } else {
    p$mat <- rbind(p$mat,simple_triplet_zero_matrix(nrow=1,ncol=NCOL(p$mat)))
  }
  p$rhs <- c(p$rhs,rhs)
  p$dir <- c(p$dir,sense)
  if( !is.null(name) ) 
    assign(name, NROW(p$mat), envir=p$ctrs) # TODO vérifier qu'elle n'existe pas déjà
  p
})

#' Set the coefficient of a variable in a contraint of a LP or MIP given their indexes
#' 
#' \code{setCoef} sets the coefficient of a variable in a contraint of
#' a linear program (LP) or mixed integer program (MIP) given its numeric
#' indexes in the problem matrix.
#' 
#' @param p an object of class \code{lpmodeler}
#' @param v index of the variable in the problem matrix (column)
#' @param c index of the constraint in the problem matrix (row)
#' @param x value of the coefficient
#' @return An object of class \code{lpmodeler}.
#' @author Cyrille Szymanski <cnszym at gmail.com>
#' @seealso TODO
#' @export
#' @examples
#' # TODO
setCoeff <- cmpfun(function(p,v,c,x)
{
  p$mat[c,v] <- x
  p
})

#' Set the coefficient of a variable in a contraint of a LP or MIP given their names
#' 
#' \code{setPoint} sets the coefficient of a variable in a contraint of
#' a linear program (LP) or mixed integer program (MIP) given their names.
#' 
#' @param p an object of class \code{lpmodeler}
#' @param v name of the variable in the problem matrix (column)
#' @param c name of the constraint in the problem matrix (row)
#' @param x value of the coefficient
#' @return An object of class \code{lpmodeler}.
#' @author Cyrille Szymanski <cnszym at gmail.com>
#' @seealso TODO
#' @export
#' @examples
#' # TODO
setPoint <- cmpfun(function(p,v,c,x)
{
  p$mat[get(c,envir=p$ctrs),get(v,envir=p$vars)] <- x
  p
})

# getCoordConstraint <- cmpfun(function(p,c) p$ctrs[[c]])
# getCoordVariable <- cmpfun(function(p,v) p$vars[[v]])

#' Solve a LP or a MIP
#' 
#' Solve a linear program (LP) or a mixed integer program (MIP).
#' Find the values of the objective function and the associated
#' variables using the specified solver.
#' 
#' @param p an object of class \code{lpmodeler}
#' @param solver name of the solver to use: Rsymphony (default)
#' @param ... other parameters passed to the solver
#' @return The object returned by the solver
#' @author Cyrille Szymanski <cnszym at gmail.com>
#' @seealso TODO
#' @export
#' @examples
#' # create and solve the following linear program:
#' # Simple mixed integer linear program.
#' # max:   3 x1 + 1 x2 + 3 x3
#' # s.t.: -1 x1 + 2 x2 +   x3 <= 4
#' #               4 x2 - 3 x3 <= 2
#' #          x1 - 3 x2 + 2 x3 <= 3
#' #          x1 >= 0 (integer)
#' #          x2 >= 0 (real)
#' #          x3 >= 0 (integer)
#' p <- newProblem()
#' p <- addVariable(p, "I", 3)
#' p <- addVariable(p, "C", 1)
#' p <- addVariable(p, "I", 3)
#' p <- addConstraint(p, "<=", 4, c(-1, 2, 1))
#' p <- addConstraint(p, "<=", 2, c(0, 4, -3))
#' p <- addConstraint(p, "<=", 3, c(1, -3, 2))
#' p <- addConstraint(p, ">=", 0, c(1, 0, 0))
#' p <- addConstraint(p, ">=", 0, c(0, 1, 0))
#' p <- addConstraint(p, ">=", 0, c(0, 0, 1))
#'
#' if(require(Rsymphony))
#'   mipSolve(p)
mipSolve <- function(p, solver=c("Rsymphony"), ...)
{
  if( solver[1]=="Rsymphony" )
    mipSolve.Rsymphony(p, ...)
  else
    stop("Unknown or unsupported solver: ", solver[1])
}
