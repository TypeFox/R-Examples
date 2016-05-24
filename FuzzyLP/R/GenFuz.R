
#' @title
#' Solves a maximization (minimization) problem having fuzzy coefficients in the constraints, the 
#' objective function and/or the technological matrix.
#'
#' @rdname GFLP
#' @description
#' The goal is to solve a linear programming problem having Trapezoidal Fuzzy Numbers 
#' as coefficients in the constraints, the objective function and/or the technological matrix. 
#' \deqn{Max\, f(x)\ or\ Min\ f(x)}{Max f(x) or Min f(x)}  
#' \deqn{s.t.:\quad \widetilde{A}x<=\widetilde{b}+(1-\beta)*\widetilde{t}}{s.t.:  Ax<=b+(1-\beta)*t}
#' This function uses different ordering functions for the objective function and for the constraints.
#'
#' @param objective A vector \eqn{(c_{1}^{f}, c_{2}^{f}, ..., c_{n}^{f})}{(c1, c2, ..., cn)} of 
#' Trapezoidal Fuzzy Numbers with the objective function coefficients 
#' \eqn{f(x)=c_{1}^{f} x_1+\ldots+c_{n}^{f} x_n}{f(x)=c1*x1+\ldots+cn*xn}. Note that any of the 
#' coefficients may also be Real Numbers.
#' @param A Technological matrix containing Trapezoidal Fuzzy Numbers and/or Real Numbers.
#' @param dir Vector of strings with the direction of the inequalities, of the same length as \code{b} and \code{t}. Each element 
#' of the vector must be one of "=", ">=", "<=", "<" or ">".
#' @param b Vector with the right hand side of the constraints. \code{b} may contain Trapezoidal Fuzzy Numbers and/or Real Numbers.
#' @param t Vector with the tolerance of each constraint. \code{t} may contain Trapezoidal Fuzzy Numbers and/or Real Numbers.
#' @param maximum \code{TRUE} to maximize the objective function, \code{FALSE} to minimize the objective function.
#' @param ordf_obj Ordering function to be used in the objective function, \code{ordf_obj} must be one of 
#' \code{"Yager1"}, \code{"Yager3"}, \code{"Adamo"} or \code{"Average"}. 
#' @param ordf_obj_param Parameters need by ordf_obj function, if it needs more than one parameter, use 
#' a named vector. See \code{\link{FOLP.ordFun}} for more information about the ordering functions parameters.
#' @param ordf_res Ordering function to be used in the constraints, \code{ordf_res} must be one of 
#' \code{"Yager1"}, \code{"Yager3"}, \code{"Adamo"} or \code{"Average"}.
#' @param ordf_res_param Parameters need by ordf_res function, if it needs more than one parameter, use 
#' a named vector. See \code{\link{FOLP.ordFun}} for more information about the ordering functions parameters.
#' @param min The lower bound of the interval to take the sample.
#' @param max The upper bound of the interval to take the sample.
#' @param step The sampling step.
#' @return \code{GFLP} returns the solutions for the sampled \eqn{\beta's} if the solver has found them. 
#' If the solver hasn't found solutions for any of the \eqn{\beta's} sampled, return NULL.
#' @references Gonzalez, A. A studing of the ranking function approach through mean values. Fuzzy Sets and Systems, 35:29-41, 1990.
#' @references Cadenas, J.M. and Verdegay, J.L. Using Fuzzy Numbers in Linear Programming. IEEE Transactions on Systems, Man, and Cybernetics-Part B: Cybernetics, vol. 27, No. 6, December 1997.
#' @references Tanaka, H., Ichihashi, H. and Asai, F. A formulation of fuzzy linear programming problems based a comparison of fuzzy numbers. Control and Cybernetics, 13:185-194, 1984.
#' @export GFLP
#' @examples 
#' ## maximize:   [1,3,4,5]*x1 + x2
#' ## s.t.:         [0,2,3,3.5]*x1 +   [0,1,1,4]*x2 <= [2,2,2,3] + (1-beta)*[1,2,2,3]
#' ##               [3,5,5,6]*x1   + [1.5,2,2,3]*x2 <= 12
#' ##               x1, x2 are non-negative real numbers
#' 
#' obj <- c(TrapezoidalFuzzyNumber(1,3,4,5), 1)
#' 
#' a11 <- TrapezoidalFuzzyNumber(0,2,2,3.5)
#' a21 <- TrapezoidalFuzzyNumber(3,5,5,6)
#' a12 <- -TrapezoidalFuzzyNumber(0,1,1,4)
#' a22 <- TrapezoidalFuzzyNumber(1.5,2,2,3)
#' A <- matrix(c(a11, a21, a12, a22), nrow = 2)
#' 
#' dir <- c("<=", "<=")
#' b<-c(TrapezoidalFuzzyNumber(2,2,2,3), 12)
#' t<-c(TrapezoidalFuzzyNumber(1,2,2,3),0);
#' max <- TRUE
#' 
#' GFLP(obj, A, dir, b, t, maximum = max, ordf_obj="Yager1", ordf_res="Yager3")
#' GFLP(obj, A, dir, b, t, maximum = max, ordf_obj="Adamo", ordf_obj_param=0.5, ordf_res="Yager3")
#' GFLP(obj, A, dir, b, t, maximum = max, "Average", ordf_obj_param=c(t=3, lambda=0.5), 
#' ordf_res="Adamo", ordf_res_param = 0.5)
#' GFLP(obj, A, dir, b, t, maximum = max, ordf_obj="Average", ordf_obj_param=c(t=3, lambda=0.8), 
#' ordf_res="Yager3", min = 0, max = 1, step = 0.2)
GFLP <- function(objective, A, dir, b, t, maximum = TRUE, 
                  ordf_obj=c("Yager1","Yager3","Adamo","Average"), ordf_obj_param=NULL, 
                  ordf_res=c("Yager1","Yager3","Adamo","Average"), ordf_res_param=NULL, 
                  min=0, max=1, step=0.25){
  old_objective<-objective
  if (any(sapply(objective,class)=="TrapezoidalFuzzyNumber")) {
    ordf_obj2=switch(EXPR=ordf_obj, "Yager1"=.Yager_1, "Yager3"=.Yager_3, "Adamo"=.Adamo, "Average"=.Average)
    
    if (ordf_obj=="Adamo"){
      objective<-sapply(objective, ordf_obj2, ordf_obj_param[1], simplify="array", USE.NAMES=FALSE)
    } else if (ordf_obj=="Average"){
      objective<-sapply(objective, ordf_obj2, ordf_obj_param[1], ordf_obj_param[2], simplify="array", USE.NAMES=FALSE)
    } else {
      objective<-sapply(objective, ordf_obj2, simplify="array", USE.NAMES=FALSE)
    }
    
    #print(objective);
  }
  
  if (any(sapply(A,class)=="TrapezoidalFuzzyNumber")) fuzmatrix=T else fuzmatrix=F;
  if (any(sapply(b,class)=="TrapezoidalFuzzyNumber")) fuzb=T else fuzb=F;
  if (any(sapply(t,class)=="TrapezoidalFuzzyNumber")) fuzt=T else fuzt=F;
  
  if(fuzmatrix || fuzb || fuzt) {
    ordf_res2=switch(EXPR=ordf_res, "Yager1"=.Yager_1, "Yager3"=.Yager_3, "Adamo"=.Adamo, "Average"=.Average)
  }
  
  if(fuzmatrix) {
    aux<-dim(A)
    
    if (ordf_res=="Adamo"){
      A<-sapply(A, ordf_res2, ordf_res_param[1], simplify="array", USE.NAMES=FALSE)
    } else if (ordf_res=="Average"){
      A<-sapply(A, ordf_res2, ordf_res_param[1], ordf_res_param[2], simplify="array", USE.NAMES=FALSE)
    } else {
      A<-sapply(A, ordf_res2, simplify="array", USE.NAMES=FALSE)
    }
    
    dim(A)<-aux
    #print(A)
  }
  
  if(fuzb) {
    if (ordf_res=="Adamo"){
      b<-sapply(b, ordf_res2, ordf_res_param[1], simplify="array", USE.NAMES=FALSE)
    } else if (ordf_res=="Average"){
      b<-sapply(b, ordf_res2, ordf_res_param[1], ordf_res_param[2], simplify="array", USE.NAMES=FALSE)
    } else {
      b<-sapply(b, ordf_res2, simplify="array", USE.NAMES=FALSE)
    }
    
    #print(b)
  }
  
  if(fuzt) {
    if (ordf_res=="Adamo"){
      t<-sapply(t, ordf_res2, ordf_res_param[1], simplify="array", USE.NAMES=FALSE)
    } else if (ordf_res=="Average"){
      t<-sapply(t, ordf_res2, ordf_res_param[1], ordf_res_param[2], simplify="array", USE.NAMES=FALSE)
    } else {
      t<-sapply(t, ordf_res2, simplify="array", USE.NAMES=FALSE)
    }
    #print(t)
  }
  #print(t)
  
  sol <- FCLP.sampledBeta(objective, A, dir, b, t, min, max, step, maximum)
  
  if (!is.null(sol)){
    sol <- sol[,1:(ncol(sol)-1)] # Deletes the last column (the objective of the aux problem)
    obj <- NULL
    for (i in 1:nrow(sol)){
      obj<-c(obj,.evalObjective(old_objective,sol[i,2:3]))
      #obj<-c(obj,7)#.evalObjective(objective,sal[i,2:3]))
    }
    
    #obj=.evalObjective(objective,sol) # Calculate the real objective
    sol <- cbind(sol, objective=obj) # Adds the real objective
    #plot(obj)
  }
  sol
}

