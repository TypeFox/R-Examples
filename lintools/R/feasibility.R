
is_contradiction <- function(A,b,neq,nleq,eps){
  
  if (nrow(A)==0||ncol(A)==0) return(logical(0))

  ieq <- seq_len(neq)
  leq <- neq + seq_len(nleq)
  lt <- neq + nleq + seq_len(nrow(A)-neq-nleq)
  # rows where ai=0
  AI <- rowSums(abs(A) > eps) == 0
  # contradiction in equality
  eqs <- AI[ieq] & abs(b[ieq]) > eps
  # contradictions in inequalities a.x <= b
  ineqs <- AI[leq] & b[leq] < -eps
  # contradictions in strict inequalities a.x < b
  sineqs <- AI[lt] & b[lt] <= 0
  
  c(eqs,  ineqs,  sineqs)
}

is_tautology <- function(A,b,neq,nleq,eps){
  if(nrow(A) == 0 |ncol(A) ==0 ) return(logical(0))
  
  ieq <- seq_len(neq)
  leq <- neq + seq_len(nleq)
  lt <- neq + nleq + seq_len(nrow(A)-neq-nleq)
  # rows where ai=0
  AI <- rowSums(abs(A) > eps) == 0
  # tautology in equality
  eqs <- AI[ieq] & abs(b[ieq]) < eps
  # tautology in inequalities a.x <= b
  ineqs <- AI[leq] & abs(b[leq]) < eps
  # tautology in strict inequalities a.x < b
  sineqs <- AI[lt] & b[lt] < -eps
  
  c(eqs,  ineqs,  sineqs)
  
  
}

are_tautologies <- function(A,b,neq,nleq,eps){
  all(is_tautology(A=A,b=b,neq=neq,nleq=nleq,eps=eps))
}

# find straigtforward contradictions of the form 0 <= b or 0 == b 
has_contradiction <- function(A,b, neq, nleq, eps){
  any(is_contradiction(A,b,neq,nleq,eps))
}


#' Check feasibility of a system of linear (in)equations
#'
#' @param A [\code{numeric}] matrix
#' @param b [\code{numeric}] vector
#' @param neq [\code{numeric}] The first \code{neq} rows in \code{A} and
#'   \code{b} are treated as linear equalities. 
#' @param nleq [\code{numeric}] The \code{nleq} rows after \code{neq} are treated as
#'   inequations of the form \code{a.x<=b}. All remaining rows are treated as strict inequations
#'   of the form \code{a.x<b}.
#' @param eps [\code{numeric}] Absolute values \code{< eps} are treated as zero.
#' @param method [\code{character}] At the moment, only the 'elimination' method is implemented.
#'
#' 
#' @examples 
#' # An infeasible system:
#' # x + y == 0
#' # x > 0
#' # y > 0
#' A <- matrix(c(1,1,1,0,0,1),byrow=TRUE,nrow=3)
#' b <- rep(0,3)
#' is_feasible(A=A,b=b,neq=1,nleq=0)
#' 
#' # A feasible system:
#' # x + y == 0
#' # x >= 0
#' # y >= 0
#' A <- matrix(c(1,1,1,0,0,1),byrow=TRUE,nrow=3)
#' b <- rep(0,3)
#' is_feasible(A=A,b=b,neq=1,nleq=2)
#' 
#' @export
is_feasible <- function(A, b, neq=nrow(A), nleq=0, eps=1e-8, method="elimination"){
  fs_elimination(A=A, b=b, neq=neq, nleq=nleq, eps=eps)
}

  ## TODO: all sorts of optimizations, including:
  # - blocking
  # - check singularity of A'A of equality section (?)
  # - figure out a good variable elimination order



fs_elimination <- function(A, b, neq, nleq, eps, H=NULL, h=0){
  # check before compact, because that also removes tautologies.
  if ( has_contradiction(A=A, b=b, neq=neq, nleq=nleq, eps=eps) ) return(FALSE)
  if ( are_tautologies(A=A, b=b, neq=neq, nleq=nleq,eps=eps) ) return(TRUE)
  L <- compact(A=A,b=b,neq=neq,nleq=nleq,eps=eps)
  # quick post-compactification check to avoid extra recursion
  if ( nrow(L$A) == 0 | ncol(L$A) == 0 ) return(TRUE)
  
  L <- eliminate(L$A, L$b, neq = L$neq, nleq=L$nleq, variable=1, H=H, h=h)
  fs_elimination(A=L$A, b=L$b, neq=L$neq,nleq=L$nleq, eps=eps, H=L$H, h=L$h) 
}






