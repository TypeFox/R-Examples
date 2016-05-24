
#' Eliminate a variable from a set of edit rules
#' 
#' Eliminating a variable amounts to deriving all (non-redundant) linear
#' (in)equations not containing that variable. Geometrically, it can be interpreted as
#' a projection of the solution space (vectors satisfying all equations) along the
#' eliminated variable's axis. 
#' 
#' @section Details:
#' For equalities Gaussian elimination is applied. If inequalities are involved,
#' Fourier-Motzkin elimination is used. In principle, FM-elimination can
#' generate a large number of redundant inequations, especially when applied
#' recursively. Redundancies can be recognized by recording how new inequations
#' have been derived from the original set. This is stored in the \code{H} matrix
#' when multiple variables are to be eliminated (Kohler, 1967).
#' 
#' 
#'
#' @param A \code{[numeric]} Matrix 
#' @param b \code{[numeric]} vector
#' @param neq [\code{numeric}] The first \code{neq} rows in \code{A} and
#'   \code{b} are treated as linear equalities. 
#' @param nleq [\code{numeric}] The \code{nleq} rows after \code{neq} are treated as
#' inequations of the form \code{a.x<=b}. All remaining rows are treated as strict inequations
#' of the form \code{a.x<b}.
#' @param variable \code{[numeric|logical|character]} Index in columns of \code{A}, representing the variable to eliminate.
#' @param eps \code{[numeric]} Coefficients with absolute value  \code{<= eps} are treated as zero.
#' @param H \code{[numeric]} (optional) Matrix indicating how linear inequalities have been derived. 
#' @param h \code{[numeric]} (optional) number indicating how many variables have been eliminated from the original system
#' using Fourier-Motzkin elimination.
#' 
#'   
#' @export
#'
#' 
#' @return A \code{list} with the folowing components
#' \itemize{
#'   \item{\code{A}: the \code{A} corresponding to the system with variables eliminated.}
#'   \item{\code{b}: the constant vector corresponding to the resulting system}
#'   \item{\code{neq}: the number of equations}
#'   \item{\code{H}: The memory matrix storing how each row was derived}
#'   \item{\code{h}: The number of variables eliminated from the original system.}
#' }
#' 
#'
#' @references
#' 
#' D.A. Kohler (1967) Projections of convex polyhedral sets, Operational Research
#' Center Report , ORC 67-29, University of California, Berkely.
#' 
#' H.P. Williams (1986) Fourier's method of linear programming and its dual. American
#' Mathematical Monthly 93, pp 681-695.
#' 
#' @examples 
#' 
#' # Example from Williams (1986)
#' A <- matrix(c(
#'    4, -5, -3,  1,
#'   -1,  1, -1,  0,
#'    1,  1,  2,  0,
#'   -1,  0,  0,  0,
#'    0, -1,  0,  0,
#'    0,  0, -1,  0),byrow=TRUE,nrow=6) 
#' b <- c(0,2,3,0,0,0)
#' L <- eliminate(A=A, b=b, neq=0, nleq=6, variable=1)
#' 
#' 
#' @export
eliminate <- function(A, b, neq=nrow(A), nleq=0, variable, H=NULL, h=0, eps=1e-8){
  check_sys(A=A, b=b, neq=neq)
 
  # TODO: use equality elimination rule.
  
  Ab <- cbind(A,b)
  if (is.character(variable)){
    var <- match(variable, colnames(A))[1]
  } else {
    var <- variable
  }
  
  ops <- rep("<",nrow(A))
  ops[seq_len(neq)] <- "=="
  ops[neq + seq_len(nleq)] <- "<="
  
  coefs <- Ab[,var]
  I <- abs(coefs) > eps
  
  eq <- I & seq_len(nrow(A)) <= neq
  
  upper <- which(!eq & coefs > 0)
  lower <- which(!eq & coefs < 0)
  
  coefs[!eq] <- abs(coefs[!eq])
  
  eq <- which(eq);
  
  # elimination possible?
  if ( (length(upper) > 0 && length(lower) > 0) ||
      (length(eq) >= 1 && (length(upper) > 0 || length(lower) > 0)) ||
      (length(eq) >= 2) ){
    h <- h+1
  } else {
    # return rows and columns where 'var' does not occur
    ii <- abs(A[,var]) < eps
    return(list(
      A = Ab[ii,-ncol(Ab), drop=FALSE]
        , b = b[ii]
        , neq = sum(which(ii)<=neq)
        , nleq = sum(which(ii)>neq & which(ii)<=neq+nleq)
        , H = H
        , h = h
    ))
  } 
  
  if ( is.null(H) ){
    H <- matrix(FALSE,nrow=nrow(A),ncol=nrow(A))
    diag(H) <- TRUE 
    colnames(H) <- rownames(A)
  }
  
  
  #normalize matrix, every row has coefficient 1 or -1 for var
  Ab[I,] <- Ab[I,] / coefs[I]
  
  # eqs and ineqs w/coef>0 ==> ineqs w/coef<0
  equpper <- c(eq, upper)
  I1 <- rep(equpper,each=length(lower))
  I2 <- rep(lower, times=length(equpper))
  ml <- Ab[I1,,drop=FALSE] + Ab[I2,,drop=FALSE]
  ol <- ifelse(ops[I1] != "<", ops[I2], ops[I1])
  dl <- H[I1,,drop=FALSE] | H[I2,,drop=FALSE]
  
  # eqs ==> ineqs w/coef>0
  I1 <- rep(eq,each=length(upper))
  I2 <- rep(upper, times=length(eq))
  mu <- Ab[I2,,drop=FALSE] - Ab[I1,,drop=FALSE]
  ou <- ops[I2]
  du <- H[I1,,drop=FALSE] | H[I2,,drop=FALSE]
  
  # eqs ==> eqs
  me <- Ab[logical(0),,drop=FALSE]
  de <- H[logical(0),,drop=FALSE]
  if ( length(eq)>1){
    me <- t(t(Ab[eq[-1],,drop=FALSE]) - Ab[eq[1],])
    de <- t(t(H[eq[-1],,drop=FALSE]) | H[eq[1],])       
  } 
  oe <- rep("==",nrow(me))
  
  Ab <- rbind(ml,mu,me,Ab[!I,,drop=FALSE])
  H <- rbind(dl,du,de,H[!I,,drop=FALSE])
  o <- c(ol,ou,oe,ops[!I])
  redundant <- rowSums(H) > h + 1 #| isObviouslyRedundant.matrix(E=m, operators=o)
  
  Ab <- Ab[!redundant,,drop=FALSE]
  H <- H[!redundant,,drop=FALSE]
  L <- normalize(
      A = Ab[,-ncol(Ab),drop=FALSE]
    , b = Ab[,ncol(Ab),drop=TRUE]
    , operators = o[!redundant]
  ) 
  list(A = L$A, b = L$b, neq=L$neq, nleq=L$nleq, H=H[L$order,,drop=FALSE], h=h)
}





