
#' @title
#' Solves a crisp linear programming problem.
#'
#' @description
#' \code{crispLP} use the classic solver (simplex) to solve a crisp linear programming problem: 
#' \deqn{Max\, f(x)\ or\ Min\ f(x)}{Max f(x) or Min f(x)}  
#' \deqn{s.t.:\quad Ax<=b}{s.t.:  Ax<=b}
#'
#' @param objective A vector \eqn{(c_1, c_2, \ldots, c_n)}{(c1, c2, \ldots, cn)} with the objective function coefficients \eqn{f(x)=c_1 x_1+\ldots+c_n x_n}{f(x)=c1*x1+\ldots+cn*xn}.
#' @param A Technological matrix of Real Numbers.
#' @param dir Vector of strings with the direction of the inequalities, of the same length as \code{b}. Each element 
#' of the vector must be one of "=", ">=", "<=", "<" or ">".
#' @param b Vector with the right hand side of the constraints.
#' @param maximum \code{TRUE} to maximize the objective function, \code{FALSE} to minimize the objective function.
#' @param verbose \code{TRUE} to show aditional screen info, \code{FALSE} to hide aditional screen info.
#' @return \code{crispLP} returns the solution if the solver has found it or NULL if not.
#' @import ROI
#' @export crispLP
#' @examples 
#' ## maximize:   3*x1 + x2
#' ## s.t.:       1.875*x1   - 1.5*x2 <= 4
#' ##              4.75*x1 + 2.125*x2 <= 14.5
#' ##               x1, x2 are non-negative real numbers
#' 
#' obj <- c(3, 1)
#' A <- matrix(c(1.875, 4.75, -1.5, 2.125), nrow = 2)
#' dir <- c("<=", "<=")
#' b <- c(4, 14.5)
#' max <- TRUE
#' 
#' crispLP(obj, A, dir, b, maximum = max, verbose = TRUE)
#' 

crispLP <- function(objective, A, dir, b, maximum = TRUE, verbose = TRUE){
  lp <- OP(objective, L_constraint(A, dir, rhs = b), maximum = maximum)
  sol <- ROI_solve(lp, solver = "glpk")
  if (verbose) print(sol$status$msg$message)
  if (sol$status$code!=0) return(NULL)
  matrix(c(sol$solution, sol$objval), ncol = length(objective)+1, byrow=TRUE, dimnames=list(NULL, c(paste("x",1:length(objective), sep=""), "objective")))
}



#' @title
#' Solves a Fuzzy Linear Programming problem with fuzzy constraints.
#'
#' @rdname FCLP.Beta
#' @description
#' The goal is to solve a linear programming problem having fuzzy constraints. 
#' \deqn{Max\, f(x)\ or\ Min\ f(x)}{Max f(x) or Min f(x)}  
#' \deqn{s.t.:\quad Ax<=b+(1-\beta)*t}{s.t.:  Ax<=b+(1-\beta)*t}
#' Where \eqn{t} means we allow not to satisfy the constraint, exceeding the bound \eqn{b} at most in \eqn{t}.
#'
#' \code{FCLP.fixedBeta} uses the classic solver (simplex) to solve the problem with a fixed value of \eqn{\beta}.
#'
#' @param objective A vector \eqn{(c_1, c_2, \ldots, c_n)}{(c1, c2, \ldots, cn)} with the objective function coefficients \eqn{f(x)=c_1 x_1+\ldots+c_n x_n}{f(x)=c1*x1+\ldots+cn*xn}.
#' @param A Technological matrix of Real Numbers.
#' @param dir Vector of strings with the direction of the inequalities, of the same length as \code{b} and \code{t}. Each element 
#' of the vector must be one of "=", ">=", "<=", "<" or ">".
#' @param b Vector with the right hand side of the constraints.
#' @param t Vector with the tolerance of each constraint.
#' @param beta The value of \eqn{\beta} to be used.
#' @param maximum \code{TRUE} to maximize the objective function, \code{FALSE} to minimize the objective function.
#' @param verbose \code{TRUE} to show aditional screen info, \code{FALSE} to hide aditional screen info.
#' @return \code{FCLP.fixedBeta} returns the solution for the given beta if the solver has found it or NULL if not.
#' @references Verdegay, J.L. Fuzzy mathematical programming. In: Fuzzy Information and Decision Processes, pages 231-237, 1982. M.M. Gupta and E.Sanchez (eds).
#' @references Delgado, M. and Herrera, F. and Verdegay, J.L. and Vila, M.A. Post-optimality analisys on the membership function of a fuzzy linear programming problem. Fuzzy Sets and Systems, 53:289-297, 1993.
#' @seealso \code{\link{FCLP.classicObjective}}, \code{\link{FCLP.fuzzyObjective}}
#' @seealso \code{\link{FCLP.fuzzyUndefinedObjective}}, \code{\link{FCLP.fuzzyUndefinedNormObjective}}
#' @export FCLP.fixedBeta
#' @examples 
#' ## maximize:   3*x1 + x2
#' ## s.t.:       1.875*x1   - 1.5*x2 <= 4 + (1-beta)*5
#' ##              4.75*x1 + 2.125*x2 <= 14.5 + (1-beta)*6
#' ##               x1, x2 are non-negative real numbers
#' 
#' obj <- c(3, 1)
#' A <- matrix(c(1.875, 4.75, -1.5, 2.125), nrow = 2)
#' dir <- c("<=", "<=")
#' b <- c(4, 14.5)
#' t <- c(5, 6)
#' valbeta <- 0.5
#' max <- TRUE
#' 
#' FCLP.fixedBeta(obj, A, dir, b, t, beta=valbeta, maximum = max, verbose = TRUE)
FCLP.fixedBeta <- function(objective, A, dir, b, t, beta=0.5, maximum = TRUE, verbose = TRUE){
  # Ax<=b+(1-beta)t , Ax>=b-(1-beta)t
  # <= -> +, >= -> -
  if (beta>1 || beta<0){
    stop("beta must be in [0,1], please use a right value", call.=FALSE)
  }
  
  out<-crispLP(objective, A, dir, b+(1-beta)*.toleranceSign(dir)*t, maximum = maximum, verbose = verbose)
  if (is.null(out)) return(NULL)
  cbind(beta=beta,out)
}



#' @title
#' Solves a Fuzzy Linear Programming problem with fuzzy constraints.
#'
#' @rdname FCLP.Beta
#' @description
#' \code{FCLP.sampledBeta} solves the problem in the same way than \code{\link{FCLP.fixedBeta}} but 
#' using several \eqn{\beta's} taking values in a sample of the \eqn{[0,1]} inteval.
#'
#' @param min The lower bound of the interval to take the sample.
#' @param max The upper bound of the interval to take the sample.
#' @param step The sampling step.
#' @return \code{FCLP.sampledBeta} returns the solutions for the sampled \eqn{\beta's} if the solver has found them. 
#' If the solver hasn't found solutions for any of the \eqn{\beta's} sampled, return NULL.
#' @export FCLP.sampledBeta
#' @examples 
#' 
#' FCLP.sampledBeta(obj, A, dir, b, t, min=0, max=1, step=0.25, maximum = max, verbose = TRUE)
FCLP.sampledBeta <- function(objective, A, dir, b, t, min=0, max=1, step=0.25, maximum = TRUE, verbose = TRUE){
  stopexec<-F
  if (min>max){
    warning("max must be greater than min, please use right values", call.=FALSE)
    stopexec<-T
  }
  if (min>1 || min<0){
    warning("min must be in [0,1], please use a right value", call.=FALSE)
    stopexec<-T
  }
  if (max>1 || max<0){
    warning("max must be in [0,1], please use a right value", call.=FALSE)
    stopexec<-T
  }
  if (step<=0){
    warning("step must be a positive value, please use right values", call.=FALSE)
    stopexec<-T
  }
  if (stopexec) stop()
  if (min==max){
    warning("To use only a value for beta you may better use FCLP.fixedFeta, the run will continue", call.=FALSE)
  }
  
  solution<-NULL
  beta<-seq(min,max,step)
  for (i in 1:length(beta)){
    sol<-FCLP.fixedBeta(objective, A, dir, b, t, beta[i], maximum = maximum, verbose = FALSE)
    if (!is.null(sol)){
      solution <- rbind(solution,sol)
    }
  }
  if (verbose && is.null(solution)) print ("Solutions not found.")
  solution
}



#' @title
#' Solves a Fuzzy Linear Programming problem with fuzzy constraints trying to assure a minimum (maximum) value 
#' of the objective function.
#'
#' @rdname FCLP.Obj
#' @description
#' The goal is to solve a linear programming problem having fuzzy constraints trying to assure a minimum (or maximum) value of the objective function.
#' \deqn{Max\, f(x)\ or\ Min\ f(x)}{Max f(x) or Min f(x)} 
#' \deqn{s.t.:\quad Ax<=b+(1-\beta)*t}{s.t.:  Ax<=b+(1-\beta)*t}
#' Where \eqn{t} means we allow not to satisfy the constraint, exceeding the bound \eqn{b} at most in \eqn{t}.
#'
#' \code{FCLP.classicObjective} solves the problem trying to assure a minimum (maximum) value \eqn{z_0}{z0} 
#' of the objective function (\eqn{f(x)>=z_0}{f(x)>=z0} in maximization problems, \eqn{f(x)<=z_0}{f(x)<=z0} in minimization 
#' problems).
#'
#' @param objective A vector \eqn{(c_1, c_2, \ldots, c_n)}{(c1, c2, \ldots, cn)} with the objective function coefficients \eqn{f(x)=c_1 x_1+\ldots+c_n x_n}{f(x)=c1*x1+\ldots+cn*xn}.
#' @param A Technological matrix of Real Numbers.
#' @param dir Vector of strings with the direction of the inequalities, of the same length as \code{b} and \code{t}. Each element 
#' of the vector must be one of "=", ">=", "<=", "<" or ">".
#' @param b Vector with the right hand side of the constraints.
#' @param t Vector with the tolerance of each constraint.
#' @param z0 The minimum (maximum in a minimization problem) value of the objective function to reach. Only 
#' used in \code{FCLP.classicObjective} and \code{FCLP.fuzzyObjective}.
#' @param maximum \code{TRUE} to maximize the objective function, \code{FALSE} to minimize the objective function.
#' @param verbose \code{TRUE} to show aditional screen info, \code{FALSE} to hide aditional screen info.
#' @return \code{FCLP.classicObjective} returns a solution reaching the given minimum (maximum) 
#' value of the objective function if the solver has found it (trying to maximize \eqn{\beta}) or NULL 
#' if not. Note that the found solution may not be the optimum for the \eqn{\beta} returned, trying \eqn{\beta} in 
#' \code{\link{FCLP.fixedBeta}} may obtain better results.
#' @references Zimmermann, H. Description and optimization of fuzzy systems. International Journal of General Systems, 2:209-215, 1976.
#' @seealso \code{\link{FCLP.fixedBeta}}, \code{\link{FCLP.sampledBeta}}
#' @export FCLP.classicObjective
#' @examples 
#' ## maximize:   3*x1 + x2 >= z0
#' ## s.t.:       1.875*x1   - 1.5*x2 <= 4 + (1-beta)*5
#' ##              4.75*x1 + 2.125*x2 <= 14.5 + (1-beta)*6
#' ##               x1, x2 are non-negative real numbers
#' 
#' obj <- c(3, 1)
#' A <- matrix(c(1.875, 4.75, -1.5, 2.125), nrow = 2)
#' dir <- c("<=", "<=")
#' b <- c(4, 14.5)
#' t <- c(5, 6)
#' max <- TRUE
#' 
#' # Problem with solution
#' FCLP.classicObjective(obj, A, dir, b, t, z0=11, maximum=max, verbose = TRUE)
#' 
#' # This problem has a bound impossible to reach
#' FCLP.classicObjective(obj, A, dir, b, t, z0=14, maximum=max, verbose = TRUE)
#' 
FCLP.classicObjective <- function(objective, A, dir, b, t, z0=0, maximum = TRUE, verbose = TRUE){
  FCLP.fuzzyObjective(objective, A, dir, b, t, z0, 0, maximum = maximum, verbose = verbose)
}


#' @title
#' Solves a Fuzzy Linear Programming problem with fuzzy constraints trying to assure a minimum (maximum) value
#' of the objective function.
#'
#' @rdname FCLP.Obj
#' @description
#' \code{FCLP.fuzzyObjective} solves the problem trying to assure a minimum (maximum) value
#' \eqn{z_0}{z0} of the objective function with tolerance \eqn{t_0}{t0} (\eqn{f(x)>=z_0-(1-\beta)t_0}{f(x)>=z0-(1-\beta)*t0} in maximization
#'  problems, \eqn{f(x)<=z_0+(1-\beta)t_0}{f(x)<=z0+(1-\beta)*t0} in minimization problems).
#'
#' @param t0 The tolerance value to the minimum (or maximum) bound for the objective function. Only 
#' used in \code{FCLP.fuzzyObjective}.
#' @return \code{FCLP.fuzzyObjective} returns a solution reaching the given minimum (maximum) 
#' value of the objective function if the solver has found it (trying to maximize \eqn{\beta}) or NULL 
#' if not. Note that the found solution may not be the optimum for the \eqn{\beta} returned, trying \eqn{\beta} in 
#' \code{\link{FCLP.fixedBeta}} may obtain better results.
#' @export FCLP.fuzzyObjective
#' @examples 
#' 
#' # This problem has a fuzzy bound impossible to reach
#' FCLP.fuzzyObjective(obj, A, dir, b, t, z0=14, t0=1, maximum=max, verbose = TRUE)
#' 
#' # This problem has a fuzzy bound reachable
#' FCLP.fuzzyObjective(obj, A, dir, b, t, z0=14, t0=2, maximum=max, verbose = TRUE)
#' 
FCLP.fuzzyObjective <- function(objective, A, dir, b, t, z0=0, t0=0, maximum = TRUE, verbose = TRUE){
  # Add objective as constraint c1x1+...cnxn> = z0-(1-beta)t0 -> (c1,...,cn)
  #1-beta=z is a new variable to minimize with constraint (1-beta)<=1 -> (0,0,...,0,1)
  A<-rbind(A,objective,0) # Add the objective as a row at the end and then a zero row
  b<-c(b,z0,1) # Add the b for the new constraints
  t<-c(t,t0,0) # Add the t for the new constraints
  
  # Add the directions of the new constraints
  if (maximum==TRUE) dir<-c(dir,">=","<=")
  else dir<-c(dir,"<=","<=")
  
  # To add a new variable (1-beta) we need a new column
  A<-cbind(A,0) # Add a zero column
  A[dim(A)[1],dim(A)[2]]<-1

  # With the new variable we need to introduce t in A
  # Ax<=b+(1-beta)t -> Ax-t(1-beta)<=b
  # Ax>=b-(1-beta)t -> Ax+t(1-beta)>=b
  A[,dim(A)[2]]<-A[,dim(A)[2]]-.toleranceSign(dir)*t
  
  # New objective is to minimize the new variable (1-beta)
  objective2<-c(objective-objective,1) # (0,0,...,0,1)
  
  sol<-crispLP(objective2, A, dir, b, maximum = FALSE, verbose = FALSE) # Solve the transformed problem (minimize 1-beta)
  if(!is.null(sol)){
    # sol has format "x1 ... xn xn+1(=1-beta) objective(=1-beta)"
    
    # Format the out like "beta x1 ... xn objective"
    # The modified last column goes to the begin 1-beta -> beta
    # After that we put x1 ... xn
    # Finally we put the last column which must be modified
    sol2<-c(beta=(1-sol[1,ncol(sol)]),sol[1,1:(ncol(sol)-2)],sol[1,ncol(sol)])
    names(sol2)[1]<-"beta"
    sol2[length(sol2)]<-sum(sol[1:(ncol(sol)-2)]*objective) # Calculate the original target
    if (verbose) cat("Bound reached, FCLP.fixedBeta with beta =",sol2[1],"may obtain better results.\n")
    return(t(as.matrix(sol2)))
  } else {
    if (maximum){
      if (verbose) print("Minimal bound not reached.")
      return(NULL)
    } else {
      if (verbose) print("Maximal bound not reached.")
      return(NULL)
    }
  }
}


#' @title
#' Solves a Fuzzy Linear Programming problem with fuzzy constraints trying to assure a minimum (maximum) 
#' value of the objective function.
#'
#' @rdname FCLP.Obj
#' @description
#' \code{FCLP.fuzzyUndefinedObjective} solves the problem trying to assure a minimum (maximum) 
#' value of the objective function with tolerance but the user doesn't fix the bound nor the 
#' tolerance. The function estimate a bound and a tolerance and call \code{\link{FCLP.fuzzyObjective}} 
#' using them.
#'
#' @return \code{FCLP.fuzzyUndefinedObjective} returns a solution reaching the estimated minimum 
#' (maximum) value of the objective function if the solver has found it (trying to maximize \eqn{\beta}) 
#' or NULL if not.
#' @references Werners, B. An interactive fuzzy programming system. Fuzzy Sets and Systems, 23:131-147, 1987.
#' @export FCLP.fuzzyUndefinedObjective
#' @examples 
#' 
#' # We want the function estimates a bound and a tolerance
#' FCLP.fuzzyUndefinedObjective(obj, A, dir, b, t, maximum=max, verbose = TRUE)
#' 
FCLP.fuzzyUndefinedObjective <- function(objective, A, dir, b, t, maximum = TRUE, verbose = TRUE){
  aux<-FCLP.sampledBeta(objective, A, dir, b, t, 0, 1, 1, maximum, verbose = FALSE) # Launch the betafixed with beta=0 and beta=1
  # Ax<=b+(1-beta)t. With beta=0, Ax<=b+t; with beta=1, Ax<=b (more restrictive)
  # Ax>=b-(1-beta)t. With beta=0, Ax>=b-t; with beta=1, Ax>=b (more restrictive)
  # The more restrictive conditions (beta=1) obtain the worse results
  if(!is.null(aux)){
    z1<-aux[[1,"objective"]] # Estimate the bound using the best results (beta=0)
    t1<-abs(z1-aux[[2,"objective"]]) # Estimate the tolerance |best results - worst results|
    if (verbose){
      print(c("Using bound = ", z1))
      print(c("Using tolerance = ", t1))
    }
    return(FCLP.fuzzyObjective(objective, A, dir, b, t, z1, t1, maximum = maximum, verbose = FALSE))
  } else
  {
    if (verbose){
      print("Solution not found.")
    }
    return(NULL)
  }
}


#' @title
#' Solves a Fuzzy Linear Programming problem with fuzzy constraints problem trying to assure a minimum (maximum) 
#' value of the objective function.
#'
#' @rdname FCLP.Obj
#' @description
#' \code{FCLP.fuzzyUndefinedNormObjective} solves the problem trying to assure a minimum (maximum) 
#' value of the objective function with tolerance but the user doesn't fix the bound nor the 
#' tolerance. The function normalize the objective, estimate a bound and a tolerance and call 
#' \code{\link{FCLP.fuzzyObjective}} using them.
#'
#' @return \code{FCLP.fuzzyUndefinedNormObjective} returns a solution reaching the estimated 
#' minimum (maximum) value of the objective function if the solver has found it (trying to 
#' maximize \eqn{\beta}) or NULL if not.
#' @references Tanaka, H. and Okuda, T. and Asai, K. On fuzzy mathematical programming. Journal of Cybernetics, 3,4:37-46, 1974.
#' @export FCLP.fuzzyUndefinedNormObjective
#' @examples 
#' 
#' # We want the function estimates a bound and a tolerance
#' FCLP.fuzzyUndefinedNormObjective(obj, A, dir, b, t, maximum=max, verbose = TRUE)
#' 
FCLP.fuzzyUndefinedNormObjective <- function(objective, A, dir, b, t, maximum = TRUE, verbose = TRUE){
  aux<-FCLP.fixedBeta(objective, A, dir, b, t, 0, maximum, verbose = FALSE) # Launch the betafixed with beta=0
  # Ax<=b+(1-beta)t. With beta=0, Ax<=b+t; with beta=1, Ax<=b (more restrictive)
  # The more restrictive conditions (beta=1) obtain the worse results
  if(!is.null(aux)){
    z1<-aux[[1,"objective"]] # Estimate the bound using the best results (beta=0)
    if (verbose) {
      print(c("Using bound = ", z1))
      print(c("Using tolerance = ", z1))
    }
    return(FCLP.fuzzyObjective(objective, A, dir, b, t, z1, z1, maximum = maximum, verbose = FALSE))
  } else
  {
    if (verbose){
      print("Solution not found.")
    }
    return(NULL)
  }
}

