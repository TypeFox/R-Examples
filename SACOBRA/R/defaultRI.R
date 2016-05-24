# defaultRI
#' 
#' Default settings for \code{\link{repairInfeasRI2}} and \code{\link{repairChootinan}}.
#' 
#' Sets suitable defaults for the repair-infeasible part of COBRA. \cr
#' With the call \code{\link{setOpts}(myRI,defaultRI())} it is possible to extend a partial list 
#' \code{myRI} to a list containing all \code{ri}-elements (the missing ones are taken from 
#' \code{defaultRI()})  
#' 
# Detail:
#' A solution \eqn{x} is said to be \strong{\eqn{\epsilon}-feasible} for constraint function \eqn{f}, if 
#'     \deqn{  f(x)+\epsilon \leq 0 }
#' \cr
#' The \strong{infeasibility} of a solution is its maximum constraint violation 
#' (0 for a feasible solution).     
#'
#'  @param repairMargin  [1e-2] repair only solutions whose infeasibility is less than this margin 
#'
#'  @return a list with the following elements:
#'    \describe{
#'      \item{RIMODE}{  [2] one out of \{0,1,2,3 \} with 0,1: deprecated older versions of RI2, 
#'          2: the recommended RI2-case, see \code{\link{repairInfeasRI2}}, 
#'          3: Chootinan's method, see \code{\link{repairChootinan}}  }
#'      \item{eps1}{  [1e-4] include all constraints not eps1-feasible into the repair mechanism  }
#'      \item{eps2}{  [1e-4] selects the solution with the shortest shift among all random 
#'          realizations which are eps2-feasible  }
#'      \item{q}{  [3.0] draw coefficients \eqn{\alpha_k} from uniform distribution \eqn{U[0,q]}  }
#'      \item{mmax}{  [1000] draw mmax random realizations  }
#'      \item{repairMargin}{  repair only solutions whose infeasibility is less than this margin. }
#'      \item{repairOnlyFresBetter}{  [FALSE] if TRUE, then repair only iterates with \cr
#'          \code{fitness < so-far-best-fitness + marFres}  }
#'      \item{marFres}{  [0.0] only relevant if \code{repairOnlyFresBetter==TRUE} }
#'    }
#'      
#'  @seealso   \code{\link{repairInfeasRI2}}, \code{\link{repairChootinan}}
#'
#'  @author Wolfgang Konen, Cologne Univeristy of Applied Sciences
#'  @export
#'
defaultRI <- function(repairMargin=1e-2) {
  ri = list(
     RIMODE=2   # 0: OLD RI, 1: RI w/o epsilon-feasibility, 2: RI2, the recommended case
                # 3: repairChootinan
    ,eps1=1e-4  # include all constraints not eps1-feasible into the repair mechanism
    ,eps2=1e-4  # selectBest() selects the solution with the shortest shift among all 
                # random realizations which are eps2-feasible 
    ,q=3.0      # draw alpha_k from uniform distribution U[0,q]
    ,mmax=1000  # draw mmax random realizations
    ,OLD=FALSE  # TRUE: activate the old repairInfeasible (before 2014-09-29) 
    ,kappa=1.2  # (only OLD) if =1.0: try to step directly to the true boundary,  
                #            if >1.0: move q bit further into the feasible region
    ,repairMargin=repairMargin    # repair only solutions whose infeasibility is less 
                # than this margin 
    ,repairOnlyFresBetter=FALSE   # if repairOnlyFresBetter=TRUE, then
                # repair only iterates with fitness < so-far-best-fitness + marFres
    ,marFres=0  # only relevant if repairOnlyFresBetter==TRUE 
  )
  return(ri);
}

# deprecated
setRI <- function(RIMODE,repairMargin=1e-2) {
  ri=defaultRI(repairMargin)
  ri$mmax=1000  # draw mmax random realizations
  
  if (RIMODE==0) ri$OLD=T
  if (RIMODE==1) ri$eps1=0
  return(ri)  
}
