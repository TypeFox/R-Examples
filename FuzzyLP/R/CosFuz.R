
#' @title
#' Solves a fuzzy objective linear programming problem using ordering functions.
#' 
#' @rdname FOLP_Ordering
#' @description
#' The goal is to solve a linear programming problem having Trapezoidal Fuzzy Numbers 
#' as coefficients in the objective function (\eqn{f(x)=c_{1}^{f} x_1+\ldots+c_{n}^{f} x_n}{f(x)=c1*x1+\ldots+cn*xn}). 
#' \deqn{Max\, f(x)\ or\ Min\ f(x)}{Max f(x) or Min f(x)}  
#' \deqn{s.t.:\quad Ax<=b}{s.t.:  Ax<=b}
#' 
#' \code{FOLP.ordFun} uses ordering functions to compare Fuzzy Numbers.
#'
#' @param objective A vector \eqn{(c_{1}^{f}, c_{2}^{f}, ..., c_{n}^{f})}{(c1, c2, ..., cn)} of 
#' Trapezoidal Fuzzy Numbers with the objective function coefficients 
#' \eqn{f(x)=c_{1}^{f} x_1+\ldots+c_{n}^{f} x_n}{f(x)=c1*x1+\ldots+cn*xn}. Note that any of the 
#' coefficients may also be Real Numbers.
#' @param A Technological matrix of Real Numbers.
#' @param dir Vector of strings with the direction of the inequalities, of the same length as \code{b}. Each element 
#' of the vector must be one of "=", ">=", "<=", "<" or ">".
#' @param b Vector with the right hand side of the constraints.
#' @param maximum \code{TRUE} to maximize the objective function, \code{FALSE} to minimize the objective function.
#' @param ordf Ordering function to be used, ordf must be one of "Yager1", "Yager3", "Adamo", "Average" or "Custom". 
#' The "Custom" option allows to use a custom linear ordering function that must be placed as FUN argument.
#' If a non linear function is used the solution may not be optimal. 
#' @param ... Additional parameters to the ordering function if needed.
#' \itemize{
#'   \item Yager1 doesn't need any parameters.
#'   \item Yager3 doesn't need any parameters.
#'   \item Adamo needs a \code{alpha} parameter which must be in the interval \code{[0,1]}.
#'   \item Average needs two parameters, \code{lambda} must be in the interval \code{[0,1]} and 
#'   \code{t} that must be greater than \code{0}.
#'   \item If Custom function needs parameters, put them here. Although not required, it is recommended 
#'   to name the parameters.
#' }
#' @param FUN Custom linear ordering function to be used if the value of ordf is "Custom". If any of the 
#' coefficients of the objective function are Real Numbers, the user must assure that the function 
#' \code{FUN} works well not only with Trapezoidal Fuzzy Numbers but also with Real Numbers.
#' @return \code{FOLP.ordFun} returns the solution if the solver has found it or NULL if not.
#' @references Gonzalez, A. A studing of the ranking function approach through mean values. Fuzzy Sets and Systems, 35:29-41, 1990.
#' @references Cadenas, J.M. and Verdegay, J.L. Using Fuzzy Numbers in Linear Programming. IEEE Transactions on Systems, Man, and Cybernetics-Part B: Cybernetics, vol. 27, No. 6, December 1997.
#' @references Tanaka, H., Ichihashi, H. and Asai, F. A formulation of fuzzy linear programming problems based a comparison of fuzzy numbers. Control and Cybernetics, 13:185-194, 1984.
#' @seealso \code{\link{FOLP.multiObj}}, \code{\link{FOLP.interv}}, \code{\link{FOLP.strat}}, \code{\link{FOLP.posib}}
#' @export FOLP.ordFun
#' @examples 
#' ## maximize:   [0,2,3]*x1 + [1,3,4,5]*x2
#' ## s.t.:         x1 + 3*x2 <= 6
#' ##               x1 +   x2 <= 4
#' ##               x1, x2 are non-negative real numbers
#' 
#' obj <- c(TrapezoidalFuzzyNumber(0,2,2,3), TrapezoidalFuzzyNumber(1,3,4,5))
#' A<-matrix(c(1, 1, 3, 1), nrow = 2)
#' dir <- c("<=", "<=")
#' b <- c(6, 4)
#' max <- TRUE
#' 
#' FOLP.ordFun(obj, A, dir, b, maximum = max, ordf="Yager1")
#' FOLP.ordFun(obj, A, dir, b, maximum = max, ordf="Yager3")
#' FOLP.ordFun(obj, A, dir, b, maximum = max, ordf="Adamo", 0.5)
#' FOLP.ordFun(obj, A, dir, b, maximum = max, ordf="Average", lambda=0.8, t=3)
#' 
#' # Define a custom linear function
#' av <- function(tfn) {mean(core(tfn))}
#' FOLP.ordFun(obj, A, dir, b, maximum = max, ordf="Custom", FUN=av)
#' 
#' # Define a custom linear function
#' avp <- function(tfn, a) {a*mean(core(tfn))}
#' FOLP.ordFun(obj, A, dir, b, maximum = max, ordf="Custom", FUN=avp, a=2)

FOLP.ordFun <- function(objective, A, dir, b, maximum = TRUE, 
                        ordf=c("Yager1","Yager3","Adamo","Average","Custom"), 
                        ..., FUN=NULL)
{
  ordf=match.arg(ordf);
  #cat("Using function... ", ordf);
  if(ordf=="Custom") {
    match.fun(FUN);
    #cat("Using custom function: ", deparse(substitute(FUN)));
  }
  indice=switch(EXPR=ordf, "Yager1"=.Yager_1, "Yager3"=.Yager_3, "Adamo"=.Adamo, "Average"=.Average, "Custom"=FUN)
  obj_ind<-sapply(objective,indice,...,simplify="array",USE.NAMES=FALSE)
  sol <- crispLP(obj_ind, A, dir, b, maximum=maximum, verbose=FALSE)

  if (!is.null(sol)){
    sol <- t(sol[,1:(ncol(sol)-1)]) # Deletes the last column (the objective of the aux problem)
    obj=.evalObjective(objective,sol) # Calculate the real objective
    sol <- cbind(sol, objective=obj) # Adds the real objective
    #plot(obj)
  }
  sol
}



#' @title
#' Solves a fuzzy objective linear programming problem using Representation Theorem.
#'
#' @rdname FOLP_Repres
#' @description
#' The goal is to solve a linear programming problem having Trapezoidal Fuzzy Numbers 
#' as coefficients in the objective function (\eqn{f(x)=c_{1}^{f} x_1+\ldots+c_{n}^{f} x_n}{f(x)=c1*x1+\ldots+cn*xn}). 
#' \deqn{Max\, f(x)\ or\ Min\ f(x)}{Max f(x) or Min f(x)}  
#' \deqn{s.t.:\quad Ax<=b}{s.t.:  Ax<=b}
#' 
#' \code{FOLP.multiObj} uses a multiobjective approach. This approach is based on each \eqn{\beta}-cut
#' of a Trapezoidal Fuzzy Number is an interval (different for each \eqn{\beta}). So the problem may be
#' considered as a Parametric Linear Problem. For a value of \eqn{\beta} fixed, the problem becomes a 
#' Multiobjective Linear Problem, this problem may be solved from different approachs, \code{FOLP.multiObj} 
#' solves it using weights, the same weight for each objective.
#'
#' @param objective A vector \eqn{(c_{1}^{f}, c_{2}^{f}, ..., c_{n}^{f})}{(c1, c2, ..., cn)} of 
#' Trapezoidal Fuzzy Numbers with the objective function coefficients 
#' \eqn{f(x)=c_{1}^{f} x_1+\ldots+c_{n}^{f} x_n}{f(x)=c1*x1+\ldots+cn*xn}. Note that any of the 
#' coefficients may also be Real Numbers.
#' @param A Technological matrix of Real Numbers.
#' @param dir Vector of strings with the direction of the inequalities, of the same length as \code{b}. Each element 
#' of the vector must be one of "=", ">=", "<=", "<" or ">".
#' @param b Vector with the right hand side of the constraints.
#' @param maximum \code{TRUE} to maximize the objective function, \code{FALSE} to minimize the objective function.
#' @param min The lower bound of the interval to take the sample.
#' @param max The upper bound of the interval to take the sample.
#' @param step The sampling step.
#' @return \code{FOLP.multiObj} returns the solutions for the sampled \eqn{\beta's} if the solver has found them. 
#' If the solver hasn't found solutions for any of the \eqn{\beta's} sampled, return NULL.
#' @references Verdegay, J.L. Fuzzy mathematical programming. In: Fuzzy Information and Decision Processes, pages 231-237, 1982. M.M. Gupta and E.Sanchez (eds).
#' @references Delgado, M. and Verdegay, J.L. and Vila, M.A. Imprecise costs in mathematical programming problems. Control and Cybernetics, 16 (2):113-121, 1987.
#' @references Bitran, G.. Linear multiple objective problems with interval coefficients. Management Science, 26(7):694-706, 1985.
#' @seealso \code{\link{FOLP.ordFun}}, \code{\link{FOLP.posib}}
#' @export FOLP.multiObj
#' @examples 
#' ## maximize:   [0,2,3]*x1 + [1,3,4,5]*x2
#' ## s.t.:         x1 + 3*x2 <= 6
#' ##               x1 +   x2 <= 4
#' ##               x1, x2 are non-negative real numbers
#' 
#' obj <- c(TrapezoidalFuzzyNumber(0,2,2,3), TrapezoidalFuzzyNumber(1,3,4,5))
#' A<-matrix(c(1, 1, 3, 1), nrow = 2)
#' dir <- c("<=", "<=")
#' b <- c(6, 4)
#' max <- TRUE
#' 
#' # Using a Multiobjective approach.
#' FOLP.multiObj(obj, A, dir, b, maximum = max, min=0, max=1, step=0.2)
#' 

FOLP.multiObj <- function(objective, A, dir, b, maximum = TRUE, min=0, max=1, step=0.25){  
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
  
  solution<-NULL
  alpha<-seq(min,max,step)
  
  for (i in 1:length(alpha)){
    MObjI<-sapply(objective,.MultiobI,alpha[i],simplify="array",USE.NAMES=FALSE)
    MObjS<-sapply(objective,.MultiobS,alpha[i],simplify="array",USE.NAMES=FALSE)
    sol<-crispLP((MObjI+MObjS)/2, A, dir, b, maximum=maximum, verbose=FALSE)
    if (!is.null(sol)){
      sol <- t(sol[,1:(ncol(sol)-1)]) # Deletes the last column (the objective of the aux problem)
      obj=.evalObjective(objective,sol) # Calculate the real objective
      sol <- cbind(alpha=alpha[i], sol, objective=obj) # Adds alpha and the real objective
      solution <- rbind(solution, sol)
      #sol <- cbind(sol, objective=obj) # Adds the real objective
      #solution <- rbind(solution,cbind(alpha=alpha[i],sol))
    }
    #if (!is.null(sol)){
    #  solution <- rbind(solution,cbind(alpha=alpha[i],sol))
    #}
  }
  solution
}


#' @title
#' Solves a fuzzy objective linear programming problem using Representation Theorem.
#'
#' @rdname FOLP_Repres
#' @description
#' \code{FOLP.interv} uses an intervalar approach. This approach is based on each \eqn{\beta}-cut
#' of a Trapezoidal Fuzzy Number is an interval (different for each \eqn{\beta}). Fixing an \eqn{\beta}, 
#' using interval arithmetic and defining an order relation for intervals is posible to compare intervals, 
#' this transforms the problem in a biobjective problem (involving the minimum and the center of intervals). 
#' Finally \code{FOLP.interv} use weights to solve the biobjective problem.
#'
#' @param w1 Weight to be used, \code{w2} is calculated as \code{w2=1-w1}. \code{w1} must 
#' be in the interval \code{[0,1]}. Only used in \code{FOLP.interv}.
#' @return \code{FOLP.interv} returns the solutions for the sampled \eqn{\beta's} if the solver has found them. 
#' If the solver hasn't found solutions for any of the \eqn{\beta's} sampled, return NULL.
#' @references Alefeld, G. and Herzberger, J. Introduction to interval computation. 1984.
#' @references Moore, R. Method and applications of interval analysis, volume 2. SIAM, 1979.
#' @export FOLP.interv
#' @examples
#' 
#' # Using a Intervalar approach.
#' FOLP.interv(obj, A, dir, b, maximum = max, w1=0.3, min=0, max=1, step=0.2)
#' 

FOLP.interv <- function(objective, A, dir, b, maximum = TRUE, w1=0.5, min=0, max=1, step=0.25){
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
  
  if (w1>1 || w1<0){
    warning("w1 must be in [0,1], execution will continue with w1=0.5", call.=FALSE)
    w1=0.5; w2=0.5;
  } else w2=1-w1;
  
  solution<-NULL
  alpha<-seq(min,max,step)
  
  for (i in 1:length(alpha)){
    obj<-sapply(objective,.Interv,alpha[i],w1,w2,simplify="array",USE.NAMES=FALSE)
    sol<-crispLP(obj, A, dir, b, maximum=maximum, verbose=FALSE)
    if (!is.null(sol)){
      sol <- t(sol[,1:(ncol(sol)-1)]) # Deletes the last column (the objective of the aux problem)
      obj=.evalObjective(objective,sol) # Calculate the real objective
      sol <- cbind(alpha=alpha[i], sol, objective=obj) # Adds alpha and the real objective
      solution <- rbind(solution, sol)
      #sol <- cbind(sol, objective=obj) # Adds the real objective
      #solution <- rbind(solution,cbind(alpha=alpha[i],sol))
    }
  }
  solution
}





#' @title
#' Solves a fuzzy objective linear programming problem using Representation Theorem.
#'
#' @rdname FOLP_Posib
#' @description
#' The goal is to solve a linear programming problem having Trapezoidal Fuzzy Numbers 
#' as coefficients in the objective function (\eqn{f(x)=c_{1}^{f} x_1+\ldots+c_{n}^{f} x_n}{f(x)=c1*x1+\ldots+cn*xn}). 
#' \deqn{Max\, f(x)\ or\ Min\ f(x)}{Max f(x) or Min f(x)}
#' \deqn{s.t.:\quad Ax<=b}{s.t.:  Ax<=b}
#' 
#' \code{FOLP.posib} uses a possibilistic approach. This approach is based on Trapezoidal Fuzzy Numbers 
#' arithmetic, so the whole objective may be considered as a Fuzzy Number itself. Defining a notion of 
#' maximum for this kind of numbers (a weighted average of the minimum and maximum of the support of 
#' the Trapezoidal number).
#'
#' @param objective A vector \eqn{(c_{1}^{f}, c_{2}^{f}, ..., c_{n}^{f})}{(c1, c2, ..., cn)} of 
#' Trapezoidal Fuzzy Numbers with the objective function coefficients 
#' \eqn{f(x)=c_{1}^{f} x_1+\ldots+c_{n}^{f} x_n}{f(x)=c1*x1+\ldots+cn*xn}. Note that any of the 
#' coefficients may also be Real Numbers.
#' @param A Technological matrix of Real Numbers.
#' @param dir Vector of strings with the direction of the inequalities, of the same length as \code{b}. Each element 
#' of the vector must be one of "=", ">=", "<=", "<" or ">".
#' @param b Vector with the right hand side of the constraints.
#' @param maximum \code{TRUE} to maximize the objective function, \code{FALSE} to minimize the objective function.
#' @param w1 Weight to be used, \code{w2} is calculated as \code{w2=1-w1}. \code{w1} must 
#' be in the interval \code{[0,1]}.
#' @return \code{FOLP.posib} returns the solution for the given weights if the solver has found it or NULL if not.
#' @references Dubois, D. and Prade, H. Operations in fuzzy numbers. International Journal of Systems Science, 9:613-626, 1978.
#' @seealso \code{\link{FOLP.ordFun}}, \code{\link{FOLP.multiObj}}, \code{\link{FOLP.interv}}, \code{\link{FOLP.strat}}
#' @export FOLP.posib
#' @examples 
#' ## maximize:   [0,2,3]*x1 + [1,3,4,5]*x2
#' ## s.t.:         x1 + 3*x2 <= 6
#' ##               x1 +   x2 <= 4
#' ##               x1, x2 are non-negative real numbers
#' 
#' obj <- c(TrapezoidalFuzzyNumber(0,2,2,3), TrapezoidalFuzzyNumber(1,3,4,5))
#' A<-matrix(c(1, 1, 3, 1), nrow = 2)
#' dir <- c("<=", "<=")
#' b <- c(6, 4)
#' max <- TRUE
#' 
#' FOLP.posib(obj, A, dir, b, maximum = max, w1=0.2)
#' 

FOLP.posib <- function(objective, A, dir, b, maximum = TRUE, w1=0.5){
  if (w1>1 || w1<0){
    warning("w1 must be in [0,1], execution will continue with w1=0.5", call.=FALSE)
    w1=0.5; w2=0.5;
  } else w2=1-w1;
  
  obj_pos<-sapply(objective,.Posibi,w1,w2,simplify="array",USE.NAMES=FALSE)
  sol <- crispLP(obj_pos, A, dir, b, maximum=maximum, verbose=FALSE)
  
  if (!is.null(sol)){
    sol <- t(sol[,1:(ncol(sol)-1)]) # Deletes the last column (the objective of the aux problem)
    obj=.evalObjective(objective,sol) # Calculate the real objective
    sol <- cbind(sol, objective=obj) # Adds the real objective
    #plot(obj)
  }
  sol
}


#' @title
#' Solves a fuzzy objective linear programming problem using Representation Theorem.
#'
#' @rdname FOLP_Repres
#' @description
#' \code{FOLP.strat} uses a stratified approach. This approach is based on that \eqn{\beta}-cuts are a 
#' sequence of nested intervals. Fixing an \eqn{\beta} two auxiliary problems are solved, the first 
#' replacing the fuzzy coefficients by the lower limits of the \eqn{\beta}-cuts, the second doing the 
#' same with the upper limits. The results of the two auxiliary problems allows to formulate a new 
#' auxiliary problem, this problem tries to maximize a parameter \eqn{\lambda}.
#' 
#' @return \code{FOLP.strat} returns the solutions and the value of \eqn{\lambda} for the sampled 
#' \eqn{\beta's} if the solver has found them. If the solver hasn't found solutions for any of the 
#' \eqn{\beta's} sampled, return NULL. A greater value of \eqn{\lambda} may be interpreted as  the 
#' obtained solution is better.
#' @references Rommelfanger, H. and Hanuscheck, R. and Wolf, J. Linear programming with fuzzy objectives. Fuzzy Sets and Systems, 29:31-48, 1989.
#' @export FOLP.strat
#' @examples
#' 
#' # Using a Stratified approach.
#' FOLP.strat(obj, A, dir, b, maximum = max, min=0, max=1, step=0.2)
#'
 
FOLP.strat <- function(objective, A, dir, b, maximum = TRUE, min=0, max=1, step=0.25){
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
  
  solution<-NULL
  alpha<-seq(min,max,step)
  nvar=ncol(A)
  nec=nrow(A)
  
  for (i in 1:length(alpha)){
    objI<-sapply(objective,.StratI,alpha[i],simplify="array",USE.NAMES=FALSE)
    solI<-crispLP(objI, A, dir, b, maximum=maximum, verbose=FALSE)    
    objS<-sapply(objective,.StratS,alpha[i],simplify="array",USE.NAMES=FALSE)
    solS<-crispLP(objS, A, dir, b, maximum=maximum, verbose=FALSE)
    
    zI1 <- solI[1,"objective"]
    zI2 <- sum(objI*solS[,1:nvar])
    zS1 <- solS[1,"objective"]
    zS2 <- sum(objS*solI[,1:nvar])
    
    #cat("zI1=",zI1," zI2=",zI2," zS1=",zS1," zS2=",zS2,"  -->  ");
    
    newobj <- c(rep(0,nvar),1)
    newA <- cbind(A, 0)
    newA <- rbind(newA, c(-objI,zI1-zI2), c(-objS,zS1-zS2))
    #newA <- rbind(newA, c(objI,-zI1+zI2), c(objS,-zS1+zS2))
    newb <- c(b,-zI2,-zS2)
    #newb <- c(b,zI2,zS2)
    if (maximum==TRUE) newdir <- c(dir,"<=","<=")
    else newdir <- c(dir,">=",">=")
    
    sol<-crispLP(newobj, newA, newdir, newb, maximum=T, verbose=F)
    
    #if (!is.null(sol)){
    #  solution <- rbind(solution, c(alpha=alpha[i],sol[,1:length(sol)-1]))
    #} else {
    #  solution <- rbind(solution, c(alpha=alpha[i],rep(NA,nvar+1)))
    #}
    if (!is.null(sol)){
      sol2 <- t(sol[,1:(ncol(sol)-1)]) # Deletes the last column (the objective of the aux problem)
      obj=.evalObjective(objective,sol2) # Calculate the real objective
      solution <- rbind(solution, c(alpha=alpha[i],sol[,1:length(sol)-1],objective=obj))
    } else {
      solution <- rbind(solution, c(alpha=alpha[i],rep(NA,nvar+2)))
    }
  }
  nam=dimnames(solution)[[2]]
  nam[nvar+2]="lambda"
  dimnames(solution)[[2]]=nam
  solution
}

